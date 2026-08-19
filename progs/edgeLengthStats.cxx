//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/low_geo/measure.hxx>
#include <src/msh_lenedg.hxx>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Metris;

namespace {

struct CustomOptions {
  std::string lengthsOutput;
  std::string objectivesOutput;
  std::vector<char*> metrisArgv;
};

struct DistributionStats {
  std::size_t count = 0;
  double minimum = 0.0;
  double percentile10 = 0.0;
  double median = 0.0;
  double mean = 0.0;
  double standardDeviation = 0.0;
  double percentile90 = 0.0;
  double maximum = 0.0;
  double log1pMean = 0.0;
  double log1pStandardDeviation = 0.0;
};

struct ElementObjectiveValues {
  int elementId = -1;
  double sizeShapeTargetVolume = 0.0;
  double sizeShapeEnergy = 0.0;
  double sizeShapeValue = 0.0;
  double stepTargetVolume = 0.0;
  double stepDistanceBarrierEnergy = 0.0;
  double stepDistanceBarrierValue = 0.0;
  double stepDistanceShapeVolumeEnergy = 0.0;
  double stepDistanceShapeVolumeValue = 0.0;
};

struct ObjectiveSummary {
  double power = 1.0;
  double energySum = 0.0;
  double targetVolume = 0.0;
  double globalEnergyRoot = 0.0;
  double globalValue = 0.0;
  DistributionStats elementEnergy;
  DistributionStats elementValue;
};

struct PointwiseMetricData {
  double trace = 0.0;
  double determinant = 0.0;
  double theta = 0.0;
};

double quantile(
    const std::vector<double>& sortedValues,
    double probability) {
  const double position =
      probability * static_cast<double>(sortedValues.size() - 1);
  const std::size_t lower = static_cast<std::size_t>(std::floor(position));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
  const double weight = position - static_cast<double>(lower);
  return sortedValues[lower] * (1.0 - weight)
       + sortedValues[upper] * weight;
}

CustomOptions extractCustomOptions(int argc, char** argv) {
  CustomOptions options;
  options.metrisArgv.push_back(argv[0]);
  for (int iarg = 1; iarg < argc; ++iarg) {
    const std::string argument(argv[iarg]);
    if (argument == "--lengths-output"
        || argument == "--objectives-output") {
      if (iarg + 1 >= argc) {
        throw std::runtime_error(argument + " requires a file name");
      }
      std::string& value = argument == "--lengths-output"
                         ? options.lengthsOutput
                         : options.objectivesOutput;
      if (!value.empty()) {
        throw std::runtime_error(argument + " may be specified only once");
      }
      value = argv[++iarg];
      continue;
    }
    options.metrisArgv.push_back(argv[iarg]);
  }
  return options;
}

DistributionStats computeDistributionStats(
    const std::vector<double>& values) {
  if (values.empty()) {
    throw std::runtime_error("cannot summarize an empty distribution");
  }

  DistributionStats stats;
  stats.minimum = std::numeric_limits<double>::infinity();
  stats.maximum = -std::numeric_limits<double>::infinity();
  double sumSquaredDifferences = 0.0;
  double log1pSumSquaredDifferences = 0.0;
  std::vector<double> sortedValues;
  sortedValues.reserve(values.size());

  for (double value : values) {
    if (!std::isfinite(value) || value < 0.0) {
      throw std::runtime_error(
          "statistics contain a negative or non-finite value");
    }

    ++stats.count;
    const double count = static_cast<double>(stats.count);
    const double delta = value - stats.mean;
    stats.mean += delta / count;
    sumSquaredDifferences += delta * (value - stats.mean);

    const double log1pValue = std::log1p(value);
    const double log1pDelta = log1pValue - stats.log1pMean;
    stats.log1pMean += log1pDelta / count;
    log1pSumSquaredDifferences +=
        log1pDelta * (log1pValue - stats.log1pMean);

    stats.minimum = std::min(stats.minimum, value);
    stats.maximum = std::max(stats.maximum, value);
    sortedValues.push_back(value);
  }

  stats.standardDeviation = std::sqrt(
      sumSquaredDifferences / static_cast<double>(stats.count));
  stats.log1pStandardDeviation = std::sqrt(
      log1pSumSquaredDifferences / static_cast<double>(stats.count));
  std::sort(sortedValues.begin(), sortedValues.end());
  stats.percentile10 = quantile(sortedValues, 0.10);
  stats.median = quantile(sortedValues, 0.50);
  stats.percentile90 = quantile(sortedValues, 0.90);
  return stats;
}

void outputDistributionJson(
    std::ostream& output,
    const DistributionStats& stats) {
  output << '{'
         << "\"count\":" << stats.count << ','
         << "\"min\":" << stats.minimum << ','
         << "\"p10\":" << stats.percentile10 << ','
         << "\"median\":" << stats.median << ','
         << "\"mean\":" << stats.mean << ','
         << "\"std\":" << stats.standardDeviation << ','
         << "\"p90\":" << stats.percentile90 << ','
         << "\"max\":" << stats.maximum << ','
         << "\"log1p_mean\":" << stats.log1pMean << ','
         << "\"log1p_std\":" << stats.log1pStandardDeviation
         << '}';
}

void outputObjectiveSummaryJson(
    std::ostream& output,
    const ObjectiveSummary& summary) {
  output << '{'
         << "\"power\":" << summary.power << ','
         << "\"energy_sum\":" << summary.energySum << ','
         << "\"target_volume\":" << summary.targetVolume << ','
         << "\"global_energy_root\":" << summary.globalEnergyRoot << ','
         << "\"global_value\":" << summary.globalValue << ','
         << "\"element_energy\":";
  outputDistributionJson(output, summary.elementEnergy);
  output << ",\"element_value\":";
  outputDistributionJson(output, summary.elementValue);
  output << '}';
}

ObjectiveSummary summarizeObjective(
    double power,
    const std::vector<double>& energies,
    const std::vector<double>& values,
    const std::vector<double>& targetVolumes) {
  if (!(power > 0.0)
      || energies.size() != values.size()
      || energies.size() != targetVolumes.size()) {
    throw std::runtime_error("invalid objective statistics inputs");
  }

  ObjectiveSummary summary;
  summary.power = power;
  for (std::size_t index = 0; index < energies.size(); ++index) {
    summary.energySum += energies[index];
    summary.targetVolume += targetVolumes[index];
  }
  if (!(summary.targetVolume > 0.0)) {
    throw std::runtime_error("objective target volume is not positive");
  }
  summary.globalEnergyRoot = std::pow(summary.energySum, 1.0 / power);
  summary.globalValue = std::pow(
      summary.energySum / summary.targetVolume, 1.0 / power);
  summary.elementEnergy = computeDistributionStats(energies);
  summary.elementValue = computeDistributionStats(values);
  return summary;
}

template <class MFT>
void outputEdgeLengthStats(
    MetrisRunner& run,
    const std::string& lengthsOutput) {
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);

  msh.cleanup();
  msh.met.setSpace(MetSpace::Exp);

  intAr2 edgeVertices;
  dblAr1 lengths;
  lenStat lengthQuality;
  getLengthEdges<MFT>(
      msh, msh.get_tdim(), -1, edgeVertices, lengths, lengthQuality,
      LenTyp::GeoSiz);

  if (lengths.get_n() == 0) {
    throw std::runtime_error("mesh has no top-dimensional edges");
  }

  // Welford's algorithm gives stable population standard deviations.
  double mean = 0.0;
  double sumSquaredDifferences = 0.0;
  double logMean = 0.0;
  double logSumSquaredDifferences = 0.0;
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  std::size_t count = 0;
  std::vector<double> sortedLengths;
  sortedLengths.reserve(lengths.get_n());

  for (double length : lengths) {
    if (!(length > 0.0)) {
      throw std::runtime_error("metric edge length is not positive");
    }

    ++count;
    const double delta = length - mean;
    mean += delta / static_cast<double>(count);
    sumSquaredDifferences += delta * (length - mean);

    const double logLength = std::log(length);
    const double logDelta = logLength - logMean;
    logMean += logDelta / static_cast<double>(count);
    logSumSquaredDifferences += logDelta * (logLength - logMean);

    minimum = std::min(minimum, length);
    maximum = std::max(maximum, length);
    sortedLengths.push_back(length);
  }

  const double standardDeviation =
      std::sqrt(sumSquaredDifferences / static_cast<double>(count));
  const double logStandardDeviation =
      std::sqrt(logSumSquaredDifferences / static_cast<double>(count));
  const double geometricMean = std::exp(logMean);
  const double geometricStandardDeviationFactor =
      std::exp(logStandardDeviation);

  std::sort(sortedLengths.begin(), sortedLengths.end());
  const double percentile10 = quantile(sortedLengths, 0.10);
  const double median = quantile(sortedLengths, 0.50);
  const double percentile90 = quantile(sortedLengths, 0.90);

  if (!lengthsOutput.empty()) {
    std::ofstream output(lengthsOutput);
    if (!output.good()) {
      throw std::runtime_error(
          "cannot open edge-length output file: " + lengthsOutput);
    }
    output << "edge_length\n" << std::setprecision(17);
    for (double length : lengths) {
      output << length << '\n';
    }
  }

  // The prefix lets callers find this record even if Metris logging is enabled.
  std::cout << std::setprecision(17)
            << "EDGE_LENGTH_STATS {"
            << "\"count\":" << count << ','
            << "\"min\":" << minimum << ','
            << "\"p10\":" << percentile10 << ','
            << "\"median\":" << median << ','
            << "\"mean\":" << mean << ','
            << "\"std\":" << standardDeviation << ','
            << "\"p90\":" << percentile90 << ','
            << "\"max\":" << maximum << ','
            << "\"log_mean\":" << logMean << ','
            << "\"log_std\":" << logStandardDeviation << ','
            << "\"geometric_mean\":" << geometricMean << ','
            << "\"geometric_std_factor\":"
            << geometricStandardDeviationFactor << ','
            << "\"fraction_unit\":" << lengthQuality.prop_unit
            << "}\n";
}

template <class MFT>
void getMetricAtQuadraturePoint(
    Mesh<MFT>& msh,
    int elementId,
    int quadraturePoint,
    const double* bary,
    const double* coordinate,
    double* metric) {
  const intAr2& elements = msh.ent2poi(2);
  if (quadraturePoint < 3) {
    const int point = elements(elementId, quadraturePoint);
    for (int component = 0; component < 3; ++component) {
      metric[component] = msh.met(point, component);
    }
  } else if constexpr(std::is_same<MFT, MetricFieldAnalytical>::value) {
    msh.met.getMetPhys(
        DifVar::None, msh.met.getSpace(), coordinate, metric, NULL);
  } else {
    msh.met.getMetBary(
        AsDeg::P1, DifVar::None, msh.met.getSpace(),
        elements[elementId], 2, bary, metric, NULL);
  }
}

double metricVolumeDensity(const double* metric) {
  const double determinant =
      metric[0] * metric[2] - metric[1] * metric[1];
  if (!(determinant > 0.0)) {
    throw std::runtime_error("target metric is not positive definite");
  }
  return std::sqrt(determinant);
}

void buildRegularJacobianTranspose(
    const double* canonicalJacobianTranspose,
    double* regularJacobianTranspose) {
  for (int row = 0; row < 2; ++row) {
    for (int column = 0; column < 2; ++column) {
      regularJacobianTranspose[2 * row + column] = 0.0;
      for (int inner = 0; inner < 2; ++inner) {
        regularJacobianTranspose[2 * row + column] +=
            Constants::invtJ_0[hana::type_c<double>][2][2 * row + inner]
            * canonicalJacobianTranspose[2 * inner + column];
      }
    }
  }
}

PointwiseMetricData evaluatePointwiseMetricData(
    const double* regularJacobianTranspose,
    const double* metric) {
  const double determinantMetric =
      metric[0] * metric[2] - metric[1] * metric[1];
  const double determinantJacobian =
      regularJacobianTranspose[0] * regularJacobianTranspose[3]
      - regularJacobianTranspose[1] * regularJacobianTranspose[2];
  if (!(determinantMetric > 0.0)
      || !(std::abs(determinantJacobian) > 0.0)) {
    throw std::runtime_error(
        "objective evaluation encountered a singular metric or element");
  }

  const double row0MetricRow0 =
      metric[0] * regularJacobianTranspose[0]
          * regularJacobianTranspose[0]
      + 2.0 * metric[1] * regularJacobianTranspose[0]
          * regularJacobianTranspose[1]
      + metric[2] * regularJacobianTranspose[1]
          * regularJacobianTranspose[1];
  const double row1MetricRow1 =
      metric[0] * regularJacobianTranspose[2]
          * regularJacobianTranspose[2]
      + 2.0 * metric[1] * regularJacobianTranspose[2]
          * regularJacobianTranspose[3]
      + metric[2] * regularJacobianTranspose[3]
          * regularJacobianTranspose[3];

  PointwiseMetricData data;
  data.trace = row0MetricRow0 + row1MetricRow1;
  // det(J^T M J) = det(J)^2 det(M). This form avoids cancellation for
  // strongly anisotropic boundary-layer metrics.
  data.determinant = determinantJacobian * determinantJacobian
                   * determinantMetric;
  data.theta = std::abs(determinantJacobian)
             * std::sqrt(determinantMetric);
  if (!(data.trace > 0.0)
      || !(data.determinant > 0.0)
      || !std::isfinite(data.trace)
      || !std::isfinite(data.determinant)
      || !std::isfinite(data.theta)) {
    throw std::runtime_error(
        "objective evaluation produced invalid metric data");
  }
  return data;
}

double regularizedDistancePower(
    double squaredDistance,
    double power,
    double regularization) {
  const double value = std::pow(
      squaredDistance + regularization * regularization,
      power / 2.0) - std::pow(regularization, power);
  return std::max(0.0, value);
}

double shiftedSizeShape(const PointwiseMetricData& data) {
  const double inverseDeterminant = 1.0 / data.determinant;
  const double sizeShape = data.trace * data.trace
                         * (1.0 + inverseDeterminant * inverseDeterminant)
                         / 8.0;
  return std::max(0.0, sizeShape - 1.0);
}

double stepDistanceIntegrand(
    const PointwiseMetricData& data,
    double power,
    double regularization,
    bool shapeVolume) {
  // Recover the two positive eigenvalue logarithms without subtracting
  // nearly equal numbers in det(A). The larger eigenvalue follows from the
  // trace and determinant; log(lambda_min)=log(det(A))-log(lambda_max).
  const double discriminantSquared = std::max(
      0.0, data.trace * data.trace - 4.0 * data.determinant);
  const double maximumEigenvalue =
      0.5 * (data.trace + std::sqrt(discriminantSquared));
  if (!(maximumEigenvalue > 0.0)) {
    throw std::runtime_error("objective matrix is not positive definite");
  }
  const double maximumLog = std::log(maximumEigenvalue);
  const double minimumLog =
      std::log(data.determinant) - maximumLog;

  double squaredDistance = 0.0;
  if (shapeVolume) {
    const double meanLog = 0.5 * (minimumLog + maximumLog);
    const double centeredMinimumLog = minimumLog - meanLog;
    const double centeredMaximumLog = maximumLog - meanLog;
    const double volumeCoordinate =
        data.determinant - 1.0 / data.determinant;
    squaredDistance = centeredMinimumLog * centeredMinimumLog
                    + centeredMaximumLog * centeredMaximumLog
                    + volumeCoordinate * volumeCoordinate / 8.0;
  } else {
    squaredDistance = minimumLog * minimumLog
                    + maximumLog * maximumLog;
  }
  return regularizedDistancePower(
      squaredDistance, power, regularization);
}

template <class MFT>
std::vector<ElementObjectiveValues> evaluateElementObjectives(
    Mesh<MFT>& msh) {
  if (msh.idim != 2 || msh.get_tdim() != 2 || msh.curdeg != 1) {
    throw std::runtime_error(
        "objective statistics currently require a 2D P1 mesh");
  }

  const int sizeShapePower = msh.param->opt_pnorm;
  const double stepDistancePower = msh.param->objective_p;
  if (sizeShapePower < 1 || stepDistancePower < 1.0) {
    throw std::runtime_error(
        "objective powers must be greater than or equal to 1");
  }

  // Statistics always evaluate the paper's positive SizeShape objective,
  // independently of the objective used to generate this mesh.
  msh.param->opt_power = 1;
  msh.param->step_distance_cavity_target_average = false;

  const intAr2& elements = msh.ent2poi(2);
  std::vector<ElementObjectiveValues> results;
  results.reserve(msh.nentt(2));

  for (int elementId = 0; elementId < msh.nentt(2); ++elementId) {
    if (isdeadent(elementId, elements)) {
      continue;
    }

    ElementObjectiveValues values;
    values.elementId = elementId;

    double barycenter[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    double centroidMetric[3] = {0.0, 0.0, 0.0};
    for (int vertex = 0; vertex < 3; ++vertex) {
      const int point = elements(elementId, vertex);
      for (int component = 0; component < 3; ++component) {
        centroidMetric[component] += msh.met(point, component) / 3.0;
      }
    }

    double physicalArea = 0.0;
    if (!isvalideltP1<2, 2>(msh, elementId, NULL, &physicalArea)
        || !(physicalArea > 0.0)) {
      throw std::runtime_error(
          "objective statistics encountered an invalid element");
    }
    values.sizeShapeTargetVolume =
        physicalArea * metricVolumeDensity(centroidMetric);
    double centroidCoordinate[2];
    double centroidCanonicalJacobianTranspose[4];
    eval2<2, 1>(
        msh.coord, elements[elementId], msh.getBasis(),
        DifVar::Bary, DifVar::None, barycenter, centroidCoordinate,
        centroidCanonicalJacobianTranspose, NULL);
    double centroidRegularJacobianTranspose[4];
    buildRegularJacobianTranspose(
        centroidCanonicalJacobianTranspose,
        centroidRegularJacobianTranspose);
    const PointwiseMetricData centroidData = evaluatePointwiseMetricData(
        centroidRegularJacobianTranspose, centroidMetric);
    const double sizeShapeValue = shiftedSizeShape(centroidData);
    values.sizeShapeEnergy = values.sizeShapeTargetVolume
                           * std::pow(sizeShapeValue, sizeShapePower);
    values.sizeShapeValue = sizeShapeValue;

    constexpr int quadraturePointCount = 4;
    for (int quadraturePoint = 0;
         quadraturePoint < quadraturePointCount;
         ++quadraturePoint) {
      double bary[3] = {0.0, 0.0, 0.0};
      if (quadraturePoint < 3) {
        bary[quadraturePoint] = 1.0;
      } else {
        for (double& coordinate : bary) {
          coordinate = 1.0 / 3.0;
        }
      }

      double coordinate[2];
      double canonicalJacobianTranspose[4];
      eval2<2, 1>(
          msh.coord, elements[elementId], msh.getBasis(),
          DifVar::Bary, DifVar::None, bary, coordinate,
          canonicalJacobianTranspose, NULL);

      double metric[3];
      getMetricAtQuadraturePoint(
          msh, elementId, quadraturePoint, bary, coordinate, metric);
      double regularJacobianTranspose[4];
      buildRegularJacobianTranspose(
          canonicalJacobianTranspose, regularJacobianTranspose);

      const PointwiseMetricData pointwiseData = evaluatePointwiseMetricData(
          regularJacobianTranspose, metric);
      values.stepTargetVolume +=
          pointwiseData.theta / quadraturePointCount;

      const double barrierIntegrand = stepDistanceIntegrand(
          pointwiseData,
          stepDistancePower,
          msh.param->step_distance_regularization,
          false);
      double barrier = 0.0;
      const double rho = std::sqrt(pointwiseData.determinant);
      if (msh.param->step_distance_barrier_beta > 0.0
          && msh.param->step_distance_barrier_rho0 > 0.0
          && rho < msh.param->step_distance_barrier_rho0) {
        const double barrierArgument = std::log(
            msh.param->step_distance_barrier_rho0 / rho);
        barrier = msh.param->step_distance_barrier_beta
                * std::pow(barrierArgument, 4.0);
      }
      values.stepDistanceBarrierEnergy +=
          (barrierIntegrand * pointwiseData.theta + barrier)
          / quadraturePointCount;

      const double shapeVolumeIntegrand = stepDistanceIntegrand(
          pointwiseData,
          stepDistancePower,
          msh.param->step_distance_regularization,
          true);
      values.stepDistanceShapeVolumeEnergy +=
          shapeVolumeIntegrand * pointwiseData.theta
          / quadraturePointCount;
    }

    if (!(values.stepTargetVolume > 0.0)) {
      throw std::runtime_error("element target volume is not positive");
    }
    values.stepDistanceBarrierValue = std::pow(
        values.stepDistanceBarrierEnergy / values.stepTargetVolume,
        1.0 / stepDistancePower);
    values.stepDistanceShapeVolumeValue = std::pow(
        values.stepDistanceShapeVolumeEnergy / values.stepTargetVolume,
        1.0 / stepDistancePower);
    results.push_back(values);
  }

  return results;
}

template <class MFT>
void outputObjectiveStats(
    MetrisRunner& run,
    const std::string& objectivesOutput) {
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  const std::vector<ElementObjectiveValues> values =
      evaluateElementObjectives(msh);
  if (values.empty()) {
    throw std::runtime_error("mesh has no top-dimensional elements");
  }

  std::vector<double> sizeShapeTargetVolumes;
  std::vector<double> sizeShapeEnergies;
  std::vector<double> sizeShapeValues;
  std::vector<double> stepTargetVolumes;
  std::vector<double> barrierEnergies;
  std::vector<double> barrierValues;
  std::vector<double> shapeVolumeEnergies;
  std::vector<double> shapeVolumeValues;
  for (const ElementObjectiveValues& element : values) {
    sizeShapeTargetVolumes.push_back(element.sizeShapeTargetVolume);
    sizeShapeEnergies.push_back(element.sizeShapeEnergy);
    sizeShapeValues.push_back(element.sizeShapeValue);
    stepTargetVolumes.push_back(element.stepTargetVolume);
    barrierEnergies.push_back(element.stepDistanceBarrierEnergy);
    barrierValues.push_back(element.stepDistanceBarrierValue);
    shapeVolumeEnergies.push_back(element.stepDistanceShapeVolumeEnergy);
    shapeVolumeValues.push_back(element.stepDistanceShapeVolumeValue);
  }

  const double sizeShapePower = msh.param->opt_pnorm;
  const double stepDistancePower = msh.param->objective_p;
  const ObjectiveSummary sizeShapeSummary = summarizeObjective(
      sizeShapePower, sizeShapeEnergies, sizeShapeValues,
      sizeShapeTargetVolumes);
  const ObjectiveSummary barrierSummary = summarizeObjective(
      stepDistancePower, barrierEnergies, barrierValues,
      stepTargetVolumes);
  const ObjectiveSummary shapeVolumeSummary = summarizeObjective(
      stepDistancePower, shapeVolumeEnergies, shapeVolumeValues,
      stepTargetVolumes);

  if (!objectivesOutput.empty()) {
    std::ofstream output(objectivesOutput);
    if (!output.good()) {
      throw std::runtime_error(
          "cannot open elemental-objective output file: "
          + objectivesOutput);
    }
    output << std::setprecision(17)
           << "element_id,sizeshape_target_volume,sizeshape_energy,"
              "sizeshape_value,step_target_volume,"
              "step_distance_barrier_energy,"
              "step_distance_barrier_value,"
              "step_distance_shape_volume_energy,"
              "step_distance_shape_volume_value\n";
    for (const ElementObjectiveValues& element : values) {
      output << element.elementId << ','
             << element.sizeShapeTargetVolume << ','
             << element.sizeShapeEnergy << ','
             << element.sizeShapeValue << ','
             << element.stepTargetVolume << ','
             << element.stepDistanceBarrierEnergy << ','
             << element.stepDistanceBarrierValue << ','
             << element.stepDistanceShapeVolumeEnergy << ','
             << element.stepDistanceShapeVolumeValue << '\n';
    }
  }

  std::cout << std::setprecision(17)
            << "OBJECTIVE_STATS {\"SizeShape\":";
  outputObjectiveSummaryJson(std::cout, sizeShapeSummary);
  std::cout << ",\"StepDistanceBarrier\":";
  outputObjectiveSummaryJson(std::cout, barrierSummary);
  std::cout << ",\"StepDistanceShapeVolume\":";
  outputObjectiveSummaryJson(std::cout, shapeVolumeSummary);
  std::cout << ",\"parameters\":{"
            << "\"sizeshape_power\":" << sizeShapePower << ','
            << "\"step_distance_power\":" << stepDistancePower << ','
            << "\"step_distance_regularization\":"
            << msh.param->step_distance_regularization << ','
            << "\"barrier_rho0\":"
            << msh.param->step_distance_barrier_rho0 << ','
            << "\"barrier_beta\":"
            << msh.param->step_distance_barrier_beta
            << "}}\n";
}

template <class MFT>
void outputAllStats(
    MetrisRunner& run,
    const CustomOptions& options) {
  outputEdgeLengthStats<MFT>(run, options.lengthsOutput);
  if (!options.objectivesOutput.empty()) {
    outputObjectiveStats<MFT>(run, options.objectivesOutput);
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 1) {
      std::cerr
          << "Usage: edgeLengthStats --in mesh.meshb --met metric.solb "
             "[--lengths-output edge_lengths.csv] "
             "[--objectives-output elemental_objectives.csv] [--verb 0]\n";
      return 2;
    }

    CustomOptions options = extractCustomOptions(argc, argv);

    MetrisRunner run(
        static_cast<int>(options.metrisArgv.size()),
        options.metrisArgv.data());
    if (run.metricFE) {
      outputAllStats<MetricFieldFE>(run, options);
    } else {
      outputAllStats<MetricFieldAnalytical>(run, options);
    }
  } catch (const std::exception& exception) {
    std::cerr << "edgeLengthStats: " << exception.what() << '\n';
    return 1;
  }

  return 0;
}

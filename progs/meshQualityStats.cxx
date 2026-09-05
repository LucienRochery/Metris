//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/io_libmeshb.hxx>
#include <src/ho_constants.hxx>
#include <src/linalg/det.hxx>
#include <src/low_eval.hxx>
#include <src/low_geo/ccoef.hxx>
#include <src/low_geo/measure.hxx>
#include <src/msh_lenedg.hxx>
#include <src/quality/low_metqua.hxx>
#include <src/quality/objective_quadrature_sample.hxx>
#include <src/quality/simplex_quadrature.hxx>
#include <src/utils/CT_loop.hxx>
#include <src/utils/aux_pp_inc.hxx>
#include <src/utils/fact_pow.hxx>

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
  std::string fieldsPrefix;
  std::vector<char*> metrisArgv;
};

struct DistributionStats {
  std::size_t count = 0;
  std::size_t finiteCount = 0;
  std::size_t positiveInfinityCount = 0;
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

struct BasicStats {
  std::size_t count = 0;
  double minimum = 0.0;
  double mean = 0.0;
  double standardDeviation = 0.0;
  double maximum = 0.0;
};

struct ElementObjectiveValues {
  int elementId = -1;
  double distortion = 0.0;
  double unit = 0.0;
  double sizeShapeTargetVolume = 0.0;
  double sizeShapeEnergy = 0.0;
  double sizeShapeValue = 0.0;
  double stepTargetVolume = 0.0;
  double stepDistanceEnergy = 0.0;
  double stepDistanceValue = 0.0;
  double stepDistanceBarrierEnergy = 0.0;
  double stepDistanceBarrierValue = 0.0;
  double stepDistanceShapeVolumeEnergy = 0.0;
  double stepDistanceShapeVolumeValue = 0.0;
  double jacobianDeterminantLowerBound = 0.0;
  double normalizedJacobianLowerBound = 0.0;
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

void outputJsonNumber(std::ostream& output, double value) {
  if (std::isfinite(value)) {
    output << value;
  } else {
    output << "null";
  }
}

CustomOptions extractCustomOptions(int argc, char** argv) {
  CustomOptions options;
  options.metrisArgv.push_back(argv[0]);
  for (int iarg = 1; iarg < argc; ++iarg) {
    const std::string argument(argv[iarg]);
    if (argument == "--lengths-output"
        || argument == "--objectives-output"
        || argument == "--fields-prefix") {
      if (iarg + 1 >= argc) {
        throw std::runtime_error(argument + " requires a file name");
      }
      std::string& value = argument == "--lengths-output"
                         ? options.lengthsOutput
                         : argument == "--objectives-output"
                         ? options.objectivesOutput
                         : options.fieldsPrefix;
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
    ++stats.count;
    if (std::isnan(value) || value < 0.0) {
      throw std::runtime_error(
          "statistics contain a negative or non-finite value: "
          + std::to_string(value));
    }
    if (std::isinf(value)) {
      ++stats.positiveInfinityCount;
      continue;
    }

    ++stats.finiteCount;
    const double count = static_cast<double>(stats.finiteCount);
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

  if (stats.finiteCount == 0) {
    return stats;
  }
  stats.standardDeviation = std::sqrt(
      sumSquaredDifferences / static_cast<double>(stats.finiteCount));
  stats.log1pStandardDeviation = std::sqrt(
      log1pSumSquaredDifferences / static_cast<double>(stats.finiteCount));
  std::sort(sortedValues.begin(), sortedValues.end());
  stats.percentile10 = quantile(sortedValues, 0.10);
  stats.median = quantile(sortedValues, 0.50);
  stats.percentile90 = quantile(sortedValues, 0.90);
  return stats;
}

BasicStats computeBasicStats(const std::vector<double>& values) {
  if (values.empty()) {
    throw std::runtime_error("cannot summarize an empty distribution");
  }
  BasicStats stats;
  stats.minimum = std::numeric_limits<double>::infinity();
  stats.maximum = -std::numeric_limits<double>::infinity();
  double sumSquaredDifferences = 0.0;
  for (double value : values) {
    if (!std::isfinite(value)) {
      throw std::runtime_error("statistics contain a non-finite value");
    }
    ++stats.count;
    const double delta = value - stats.mean;
    stats.mean += delta / static_cast<double>(stats.count);
    sumSquaredDifferences += delta * (value - stats.mean);
    stats.minimum = std::min(stats.minimum, value);
    stats.maximum = std::max(stats.maximum, value);
  }
  stats.standardDeviation = std::sqrt(
      sumSquaredDifferences / static_cast<double>(stats.count));
  return stats;
}

void validateNonnegativeDistribution(
    const std::string& name,
    const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (std::isnan(values[index]) || values[index] < 0.0) {
      throw std::runtime_error(
          name + " contains an invalid value at entry "
          + std::to_string(index) + ": " + std::to_string(values[index]));
    }
  }
}

void outputDistributionJson(
    std::ostream& output,
    const DistributionStats& stats) {
  output << '{'
         << "\"count\":" << stats.count << ','
         << "\"finite_count\":" << stats.finiteCount << ','
         << "\"positive_infinity_count\":"
         << stats.positiveInfinityCount << ',';
  if (stats.finiteCount == 0) {
    output << "\"min\":null,\"p10\":null,\"median\":null,"
              "\"mean\":null,\"std\":null,\"p90\":null,"
              "\"max\":null,\"log1p_mean\":null,"
              "\"log1p_std\":null}";
    return;
  }
  output << "\"min\":";
  outputJsonNumber(output, stats.minimum);
  output << ",\"p10\":";
  outputJsonNumber(output, stats.percentile10);
  output << ",\"median\":";
  outputJsonNumber(output, stats.median);
  output << ",\"mean\":";
  outputJsonNumber(output, stats.mean);
  output << ",\"std\":";
  outputJsonNumber(output, stats.standardDeviation);
  output << ",\"p90\":";
  outputJsonNumber(output, stats.percentile90);
  output << ",\"max\":";
  outputJsonNumber(output, stats.maximum);
  output << ",\"log1p_mean\":";
  outputJsonNumber(output, stats.log1pMean);
  output << ",\"log1p_std\":";
  outputJsonNumber(output, stats.log1pStandardDeviation);
  output << '}';
}

void outputBasicStatsJson(std::ostream& output, const BasicStats& stats) {
  output << '{'
         << "\"count\":" << stats.count << ','
         << "\"min\":" << stats.minimum << ','
         << "\"mean\":" << stats.mean << ','
         << "\"std\":" << stats.standardDeviation << ','
         << "\"max\":" << stats.maximum
         << '}';
}

void outputObjectiveSummaryJson(
    std::ostream& output,
    const ObjectiveSummary& summary) {
  output << '{'
         << "\"power\":";
  outputJsonNumber(output, summary.power);
  output << ",\"energy_sum\":";
  outputJsonNumber(output, summary.energySum);
  output << ",\"target_volume\":";
  outputJsonNumber(output, summary.targetVolume);
  output << ",\"global_energy_root\":";
  outputJsonNumber(output, summary.globalEnergyRoot);
  output << ",\"global_value\":";
  outputJsonNumber(output, summary.globalValue);
  output << ",\"element_energy\":";
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

template <class MFT, int gdim, int tdim, int degree>
double integrateTargetMeasure(Mesh<MFT>& msh, const int* element) {
  const SimplexQuadratureView<tdim> quadrature =
      get_objective_quadrature<tdim>(msh.param->objective_quadrature_order);
  double measure = 0.0;
  for (int quadraturePoint = 0;
       quadraturePoint < quadrature.size();
       ++quadraturePoint) {
    const auto sample = prepare_objective_quadrature_sample<
        MFT, gdim, tdim, degree>(
            msh, AsDeg::Pk, element, quadrature[quadraturePoint],
            standard_objective_quadrature_theta_mode());
    if (!sample.theta_is_valid) {
      throw std::runtime_error(
          "target-measure quadrature encountered an invalid element");
    }
    measure += sample.quadrature_weight * sample.theta;
  }
  return measure;
}

template <int gdim, int degree>
void getJacobianLowerBounds(
    const MeshBase& msh,
    int elementId,
    double* determinantLowerBound,
    double* normalizedLowerBound) {
  constexpr int jacobianDegree = gdim * (degree - 1);
  constexpr int coefficientCount = getnnode(gdim, jacobianDegree);
  double coefficients[coefficientCount];
  getccoef<gdim, gdim, degree>(msh, elementId, NULL, coefficients);

  *determinantLowerBound = coefficients[0];
  for (int coefficient = 1;
       coefficient < coefficientCount;
       ++coefficient) {
    *determinantLowerBound = std::min(
        *determinantLowerBound, coefficients[coefficient]);
  }

  bool invalidP1 = false;
  const double p1Measure = getmeasentP1<gdim, gdim>(
      msh, elementId, NULL, &invalidP1);
  const double normalization = std::abs(p1Measure) * ifact<gdim>();
  if (!(normalization > 0.0) || !std::isfinite(normalization)) {
    throw std::runtime_error(
        "Jacobian statistics encountered a degenerate corner element");
  }
  *normalizedLowerBound = *determinantLowerBound / normalization;
}

template <int gdim, int geometryDegree, int fieldDegree>
void writeJacobianDeterminantDGField(
    const MeshBase& msh,
    const std::string& fieldsPrefix) {
  constexpr int fieldNodeCount = getnnode(gdim, fieldDegree);
  constexpr auto fieldOrdering = ORDELT(gdim);
  const intAr2& elements = msh.ent2poi(gdim);
  dblAr1 field(msh.nentt(gdim) * fieldNodeCount);

  for (int elementId = 0;
       elementId < msh.nentt(gdim);
       ++elementId) {
    if (isdeadent(elementId, elements)) {
      for (int fieldNode = 0;
           fieldNode < fieldNodeCount;
           ++fieldNode) {
        field[elementId * fieldNodeCount + fieldNode] = 0.0;
      }
      continue;
    }

    for (int fieldNode = 0;
         fieldNode < fieldNodeCount;
         ++fieldNode) {
      double barycentricCoordinates[gdim + 1];
      for (int coordinate = 0;
           coordinate < gdim + 1;
           ++coordinate) {
        barycentricCoordinates[coordinate] =
            fieldOrdering[fieldDegree][fieldNode][coordinate]
            / static_cast<double>(fieldDegree);
      }

      double physicalCoordinates[gdim];
      double jacobianTranspose[gdim * gdim];
      if constexpr (gdim == 2) {
        eval2<gdim, geometryDegree>(
            msh.coord, elements[elementId], msh.getBasis(),
            DifVar::Bary, DifVar::None, barycentricCoordinates,
            physicalCoordinates, jacobianTranspose, NULL);
      } else {
        eval3<gdim, geometryDegree>(
            msh.coord, elements[elementId], msh.getBasis(),
            DifVar::Bary, DifVar::None, barycentricCoordinates,
            physicalCoordinates, jacobianTranspose, NULL);
      }
      field[elementId * fieldNodeCount + fieldNode] =
          detmat<gdim>(jacobianTranspose);
    }
  }

  writeField(
      fieldsPrefix + ".jacobian_determinant.solb",
      msh, SolTyp::DG, field, 1, fieldDegree);
}

template <class MFT, int gdim, int degree>
std::vector<ElementObjectiveValues> evaluateElementObjectivesImpl(
    Mesh<MFT>& msh) {
  const int tdim = msh.get_tdim();
  if (tdim != gdim) {
    throw std::runtime_error(
        "objective and Jacobian statistics require a full-dimensional mesh");
  }

  const double objectivePower = msh.param->objective_p;
  if (objectivePower < 1.0) {
    throw std::runtime_error(
        "the objective power must be greater than or equal to 1");
  }

  const bool originalShapeVolume = msh.param->step_distance_shape_volume;
  const bool originalTargetAverage =
      msh.param->step_distance_cavity_target_average;
  const double originalBarrierBeta =
      msh.param->step_distance_barrier_beta;
  msh.param->step_distance_cavity_target_average = false;

  const intAr2& elements = msh.ent2poi(gdim);
  std::vector<ElementObjectiveValues> results;
  results.reserve(msh.nentt(gdim));

  for (int elementId = 0;
       elementId < msh.nentt(gdim);
       ++elementId) {
    if (isdeadent(elementId, elements)) {
      continue;
    }

    ElementObjectiveValues values;
    values.elementId = elementId;
    values.distortion = metqua<MFT, gdim, gdim, QuaFun::Distortion>(
        msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);
    values.unit = metqua<MFT, gdim, gdim, QuaFun::Unit>(
        msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);

    const double targetMeasure = integrateTargetMeasure<
        MFT, gdim, gdim, degree>(msh, elements[elementId]);
    if (!(targetMeasure > 0.0) || !std::isfinite(targetMeasure)) {
      throw std::runtime_error("element target measure is not positive");
    }
    values.sizeShapeTargetVolume = targetMeasure;
    values.stepTargetVolume = targetMeasure;

    msh.param->step_distance_shape_volume = false;
    values.sizeShapeEnergy =
        metqua<MFT, gdim, gdim, QuaFun::SizeShape>(
            msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);

    msh.param->step_distance_barrier_beta = 0.0;
    values.stepDistanceEnergy =
        metqua<MFT, gdim, gdim, QuaFun::StepDistance>(
            msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);

    msh.param->step_distance_barrier_beta = originalBarrierBeta;
    values.stepDistanceBarrierEnergy =
        metqua<MFT, gdim, gdim, QuaFun::StepDistance>(
            msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);

    msh.param->step_distance_shape_volume = true;
    values.stepDistanceShapeVolumeEnergy =
        metqua<MFT, gdim, gdim, QuaFun::StepDistance>(
            msh, AsDeg::Pk, AsDeg::Pk, elementId, 1.0);

    values.sizeShapeValue = std::pow(
        values.sizeShapeEnergy / targetMeasure, 1.0 / objectivePower);
    values.stepDistanceValue = std::pow(
        values.stepDistanceEnergy / targetMeasure,
        1.0 / objectivePower);
    values.stepDistanceBarrierValue = std::pow(
        values.stepDistanceBarrierEnergy / targetMeasure,
        1.0 / objectivePower);
    values.stepDistanceShapeVolumeValue = std::pow(
        values.stepDistanceShapeVolumeEnergy / targetMeasure,
        1.0 / objectivePower);
    getJacobianLowerBounds<gdim, degree>(
        msh, elementId,
        &values.jacobianDeterminantLowerBound,
        &values.normalizedJacobianLowerBound);
    results.push_back(values);
  }

  msh.param->step_distance_shape_volume = originalShapeVolume;
  msh.param->step_distance_cavity_target_average = originalTargetAverage;
  msh.param->step_distance_barrier_beta = originalBarrierBeta;
  return results;
}

template <class MFT>
std::vector<ElementObjectiveValues> evaluateElementObjectives(
    Mesh<MFT>& msh) {
  std::vector<ElementObjectiveValues> results;
  bool dispatched = false;
  CT_FOR0_INC(1, METRIS_MAX_DEG, degree) {
    if (degree == msh.curdeg) {
      if (msh.idim == 2) {
        results = evaluateElementObjectivesImpl<MFT, 2, degree>(msh);
      } else if (msh.idim == 3) {
        results = evaluateElementObjectivesImpl<MFT, 3, degree>(msh);
      }
      dispatched = true;
    }
  } CT_FOR1(degree);
  if (!dispatched) {
    throw std::runtime_error("unsupported mesh geometry degree");
  }
  return results;
}

template <class MFT>
void outputObjectiveStats(
    MetrisRunner& run,
    const std::string& objectivesOutput,
    const std::string& fieldsPrefix) {
  Mesh<MFT>& msh = static_cast<Mesh<MFT>&>(*run.msh_g);
  const std::vector<ElementObjectiveValues> values =
      evaluateElementObjectives(msh);
  if (values.empty()) {
    throw std::runtime_error("mesh has no top-dimensional elements");
  }

  std::vector<double> sizeShapeTargetVolumes;
  std::vector<double> distortionValues;
  std::vector<double> unitValues;
  std::vector<double> sizeShapeEnergies;
  std::vector<double> sizeShapeValues;
  std::vector<double> stepTargetVolumes;
  std::vector<double> stepEnergies;
  std::vector<double> stepValues;
  std::vector<double> barrierEnergies;
  std::vector<double> barrierValues;
  std::vector<double> shapeVolumeEnergies;
  std::vector<double> shapeVolumeValues;
  std::vector<double> jacobianLowerBounds;
  std::vector<double> normalizedJacobianLowerBounds;
  for (const ElementObjectiveValues& element : values) {
    distortionValues.push_back(element.distortion);
    unitValues.push_back(element.unit);
    sizeShapeTargetVolumes.push_back(element.sizeShapeTargetVolume);
    sizeShapeEnergies.push_back(element.sizeShapeEnergy);
    sizeShapeValues.push_back(element.sizeShapeValue);
    stepTargetVolumes.push_back(element.stepTargetVolume);
    stepEnergies.push_back(element.stepDistanceEnergy);
    stepValues.push_back(element.stepDistanceValue);
    barrierEnergies.push_back(element.stepDistanceBarrierEnergy);
    barrierValues.push_back(element.stepDistanceBarrierValue);
    shapeVolumeEnergies.push_back(element.stepDistanceShapeVolumeEnergy);
    shapeVolumeValues.push_back(element.stepDistanceShapeVolumeValue);
    jacobianLowerBounds.push_back(element.jacobianDeterminantLowerBound);
    normalizedJacobianLowerBounds.push_back(
        element.normalizedJacobianLowerBound);
  }

  const double sizeShapePower = msh.param->objective_p;
  const double stepDistancePower = msh.param->objective_p;
  validateNonnegativeDistribution("Distortion", distortionValues);
  validateNonnegativeDistribution("Unit", unitValues);
  validateNonnegativeDistribution("SizeShape energy", sizeShapeEnergies);
  validateNonnegativeDistribution("SizeShape value", sizeShapeValues);
  validateNonnegativeDistribution("StepDistance energy", stepEnergies);
  validateNonnegativeDistribution("StepDistance value", stepValues);
  validateNonnegativeDistribution(
      "StepDistanceBarrier energy", barrierEnergies);
  validateNonnegativeDistribution(
      "StepDistanceBarrier value", barrierValues);
  validateNonnegativeDistribution(
      "StepDistanceShapeVolume energy", shapeVolumeEnergies);
  validateNonnegativeDistribution(
      "StepDistanceShapeVolume value", shapeVolumeValues);
  const DistributionStats distortionSummary =
      computeDistributionStats(distortionValues);
  const DistributionStats unitSummary = computeDistributionStats(unitValues);
  const ObjectiveSummary sizeShapeSummary = summarizeObjective(
      sizeShapePower, sizeShapeEnergies, sizeShapeValues,
      sizeShapeTargetVolumes);
  const ObjectiveSummary stepSummary = summarizeObjective(
      stepDistancePower, stepEnergies, stepValues,
      stepTargetVolumes);
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
           << "element_id,distortion,unit,sizeshape_target_volume,"
              "sizeshape_energy,"
              "sizeshape_value,step_target_volume,"
              "step_distance_energy,"
              "step_distance_value,"
              "step_distance_barrier_energy,"
              "step_distance_barrier_value,"
              "step_distance_shape_volume_energy,"
              "step_distance_shape_volume_value,"
              "jacobian_determinant_lower_bound,"
              "normalized_jacobian_lower_bound\n";
    for (const ElementObjectiveValues& element : values) {
      output << element.elementId << ','
             << element.distortion << ','
             << element.unit << ','
             << element.sizeShapeTargetVolume << ','
             << element.sizeShapeEnergy << ','
             << element.sizeShapeValue << ','
             << element.stepTargetVolume << ','
             << element.stepDistanceEnergy << ','
             << element.stepDistanceValue << ','
             << element.stepDistanceBarrierEnergy << ','
             << element.stepDistanceBarrierValue << ','
             << element.stepDistanceShapeVolumeEnergy << ','
             << element.stepDistanceShapeVolumeValue << ','
             << element.jacobianDeterminantLowerBound << ','
             << element.normalizedJacobianLowerBound << '\n';
    }
  }

  if (!fieldsPrefix.empty()) {
    auto writeElementField = [&](const std::string& suffix, auto getter) {
      dblAr1 field(msh.nentt(msh.get_tdim()));
      for (int elementId = 0;
           elementId < msh.nentt(msh.get_tdim());
           ++elementId) {
        field[elementId] = 0.0;
      }
      for (const ElementObjectiveValues& element : values) {
        field[element.elementId] = getter(element);
      }
      writeField(
          fieldsPrefix + "." + suffix + ".solb",
          msh, SolTyp::P0Elt, field);
    };
    writeElementField("distortion", [](const auto& element) {
      return element.distortion;
    });
    writeElementField("unit", [](const auto& element) {
      return element.unit;
    });
    writeElementField("sizeshape_energy", [](const auto& element) {
      return element.sizeShapeEnergy;
    });
    writeElementField("sizeshape_value", [](const auto& element) {
      return element.sizeShapeValue;
    });
    writeElementField("step_distance_energy", [](const auto& element) {
      return element.stepDistanceEnergy;
    });
    writeElementField("step_distance_value", [](const auto& element) {
      return element.stepDistanceValue;
    });
    writeElementField("step_distance_barrier_energy", [](const auto& element) {
      return element.stepDistanceBarrierEnergy;
    });
    writeElementField("step_distance_barrier_value", [](const auto& element) {
      return element.stepDistanceBarrierValue;
    });
    writeElementField(
        "step_distance_shape_volume_energy", [](const auto& element) {
          return element.stepDistanceShapeVolumeEnergy;
        });
    writeElementField(
        "step_distance_shape_volume_value", [](const auto& element) {
          return element.stepDistanceShapeVolumeValue;
        });
    writeElementField(
        "jacobian_determinant_lower_bound", [](const auto& element) {
          return element.jacobianDeterminantLowerBound;
        });
    writeElementField(
        "normalized_jacobian_lower_bound", [](const auto& element) {
          return element.normalizedJacobianLowerBound;
        });
    if (msh.curdeg == 1 && msh.idim == 2) {
      writeJacobianDeterminantDGField<2, 1, 1>(msh, fieldsPrefix);
    } else if (msh.curdeg == 1 && msh.idim == 3) {
      writeJacobianDeterminantDGField<3, 1, 1>(msh, fieldsPrefix);
    } else if (msh.curdeg == 2 && msh.idim == 2) {
      writeJacobianDeterminantDGField<2, 2, 2>(msh, fieldsPrefix);
    } else if (msh.curdeg == 2 && msh.idim == 3) {
      writeJacobianDeterminantDGField<3, 2, 3>(msh, fieldsPrefix);
    } else {
      throw std::runtime_error(
          "Jacobian DG output currently supports P1/P2 geometry");
    }
    std::cout << "SOLUTION_FIELDS {\"prefix\":\""
              << fieldsPrefix << "\",\"count\":13}\n";
  }

  std::cout << std::setprecision(17)
            << "OBJECTIVE_STATS {\"Distortion\":";
  outputDistributionJson(std::cout, distortionSummary);
  std::cout << ",\"Unit\":";
  outputDistributionJson(std::cout, unitSummary);
  std::cout << ",\"SizeShape\":";
  outputObjectiveSummaryJson(std::cout, sizeShapeSummary);
  std::cout << ",\"StepDistance\":";
  outputObjectiveSummaryJson(std::cout, stepSummary);
  std::cout << ",\"StepDistanceBarrier\":";
  outputObjectiveSummaryJson(std::cout, barrierSummary);
  std::cout << ",\"StepDistanceShapeVolume\":";
  outputObjectiveSummaryJson(std::cout, shapeVolumeSummary);
  std::cout << ",\"JacobianDeterminantLowerBound\":";
  outputBasicStatsJson(std::cout, computeBasicStats(jacobianLowerBounds));
  std::cout << ",\"NormalizedJacobianLowerBound\":";
  outputBasicStatsJson(
      std::cout, computeBasicStats(normalizedJacobianLowerBounds));
  std::cout << ",\"parameters\":{"
            << "\"sizeshape_power\":" << sizeShapePower << ','
            << "\"step_distance_power\":" << stepDistancePower << ','
            << "\"classical_pnorm\":" << msh.param->opt_pnorm << ','
            << "\"geometry_degree\":" << msh.curdeg << ','
            << "\"objective_quadrature_order\":"
            << msh.param->objective_quadrature_order << ','
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
  outputObjectiveStats<MFT>(
      run, options.objectivesOutput, options.fieldsPrefix);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 1) {
      std::cerr
          << "Usage: meshQualityStats --in mesh.meshb --met metric.solb "
             "[--lengths-output edge_lengths.csv] "
             "[--objectives-output elemental_objectives.csv] "
             "[--fields-prefix mesh_quality] [--verb 0]\n";
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
    std::cerr << "meshQualityStats: " << exception.what() << '\n';
    return 1;
  }

  return 0;
}

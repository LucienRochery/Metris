//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
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

std::string extractCustomOption(
    int argc,
    char** argv,
    const std::string& option,
    std::vector<char*>& metrisArgv) {
  std::string value;
  metrisArgv.clear();
  metrisArgv.push_back(argv[0]);

  for (int iarg = 1; iarg < argc; ++iarg) {
    if (argv[iarg] == option) {
      if (iarg + 1 >= argc) {
        throw std::runtime_error(option + " requires a file name");
      }
      value = argv[++iarg];
      continue;
    }
    metrisArgv.push_back(argv[iarg]);
  }
  return value;
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

  // Welford's algorithm gives a stable population standard deviation.
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

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 1) {
      std::cerr
          << "Usage: edgeLengthStats --in mesh.meshb --met metric.solb "
             "[--lengths-output edge_lengths.csv] [--verb 0]\n";
      return 2;
    }

    std::vector<char*> metrisArgv;
    const std::string lengthsOutput =
        extractCustomOption(
            argc, argv, "--lengths-output", metrisArgv);

    MetrisRunner run(
        static_cast<int>(metrisArgv.size()), metrisArgv.data());
    if (run.metricFE) {
      outputEdgeLengthStats<MetricFieldFE>(run, lengthsOutput);
    } else {
      outputEdgeLengthStats<MetricFieldAnalytical>(run, lengthsOutput);
    }
  } catch (const std::exception& exception) {
    std::cerr << "edgeLengthStats: " << exception.what() << '\n';
    return 1;
  }

  return 0;
}

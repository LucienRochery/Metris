// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/io_libmeshb.hxx>
#include <src/low_geo/validity.hxx>
#include <src/quality/low_metqua_d.hxx>
#include <src/quality/msh_metqua.hxx>
#include <src/smoothing/msh_smooball.hxx>

#include <nlopt.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Metris;

namespace
{

enum class VariableKind
{
  PhysicalX,
  PhysicalY,
  CurveParameter
};

enum class GradientMode
{
  FrozenAnalytical,
  FiniteDifference
};

struct GeometryVariable
{
  VariableKind kind;
  int point;
  int boundaryRecord = -1;
  int cadEdgeReference = -1;
  ego cadEdge = nullptr;
};

struct CustomOptions
{
  std::filesystem::path elevatedOutput;
  std::filesystem::path historyOutput;
  QuaFun objective = QuaFun::SizeShape;
  std::string objectiveName = "SizeShape";
  GradientMode gradientMode = GradientMode::FiniteDifference;
  std::string gradientName = "finite-difference";
  int maximumEvaluations = 100;
  std::vector<char*> metrisArgv;
};

CustomOptions extractCustomOptions(int argc, char** argv)
{
  CustomOptions options;
  options.metrisArgv.push_back(argv[0]);
  for(int argument = 1; argument < argc; argument++)
  {
    const std::string name(argv[argument]);
    if(name == "--elevated-output" || name == "--history-output"
       || name == "--smoothing-objective" || name == "--gradient-mode"
       || name == "--maximum-evaluations")
    {
      if(argument + 1 >= argc)
        throw std::runtime_error(name + " requires one value");
      const std::string value(argv[++argument]);
      if(name == "--elevated-output") options.elevatedOutput = value;
      else if(name == "--history-output") options.historyOutput = value;
      else if(name == "--maximum-evaluations")
        options.maximumEvaluations = std::stoi(value);
      else if(name == "--smoothing-objective")
      {
        if(value == "sizeshape")
        {
          options.objective = QuaFun::SizeShape;
          options.objectiveName = "SizeShape";
        }
        else if(value == "step-distance-sv")
        {
          options.objective = QuaFun::StepDistance;
          options.objectiveName = "StepDistanceSV";
        }
        else
          throw std::runtime_error(
              "--smoothing-objective must be sizeshape or step-distance-sv");
      }
      else
      {
        if(value == "frozen")
        {
          options.gradientMode = GradientMode::FrozenAnalytical;
          options.gradientName = "frozen";
        }
        else if(value == "finite-difference")
        {
          options.gradientMode = GradientMode::FiniteDifference;
          options.gradientName = "finite-difference";
        }
        else
          throw std::runtime_error(
              "--gradient-mode must be frozen or finite-difference");
      }
      continue;
    }
    options.metrisArgv.push_back(argv[argument]);
  }
  if(options.elevatedOutput.empty())
    throw std::runtime_error("--elevated-output is required");
  if(options.historyOutput.empty())
    throw std::runtime_error("--history-output is required");
  if(options.maximumEvaluations < 1)
    throw std::runtime_error("--maximum-evaluations must be positive");
  return options;
}

struct GlobalSmoothingProblem
{
  Mesh<MetricFieldAnalytical>& mesh;
  QuaFun objective;
  GradientMode gradientMode;
  std::vector<GeometryVariable> variables;
  std::vector<double> initial;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> best;
  double bestValue = std::numeric_limits<double>::infinity();
  std::ofstream& history;
  std::chrono::steady_clock::time_point startTime;
  int callbackCount = 0;
  int validValueEvaluations = 0;
  int invalidValueEvaluations = 0;

  void updateAnalyticalMetricAtNodes()
  {
    mesh.met.setSpace(MetSpace::Exp);
    for(int point = 0; point < mesh.npoin; point++)
    {
      if(mesh.isdeadpoint(point)) continue;
      mesh.met.getMetPhys(
          DifVar::None,MetSpace::Exp,mesh.coord[point],mesh.met[point],nullptr);
    }
  }

  void apply(const std::vector<double>& x)
  {
    METRIS_ENFORCE(x.size() == variables.size());
    for(std::size_t index = 0; index < variables.size(); index++)
    {
      const GeometryVariable& variable = variables[index];
      if(variable.kind == VariableKind::PhysicalX)
        mesh.coord(variable.point,0) = x[index];
      else if(variable.kind == VariableKind::PhysicalY)
        mesh.coord(variable.point,1) = x[index];
      else
      {
        double evaluation[18];
        const double parameter[2] = {x[index],0.};
        const int status = EG_evaluate(variable.cadEdge,parameter,evaluation);
        if(status != EGADS_SUCCESS)
          throw std::runtime_error("EG_evaluate failed in global smoothing");
        mesh.bpo2rbi(variable.boundaryRecord,0) = x[index];
        mesh.coord(variable.point,0) = evaluation[0];
        mesh.coord(variable.point,1) = evaluation[1];
      }
    }
    updateAnalyticalMetricAtNodes();
  }

  bool geometryIsValid() const
  {
    for(int triangle = 0; triangle < mesh.nface; triangle++)
    {
      if(isdeadent(triangle,mesh.fac2poi)) continue;
      if(!classify_element_validity<2,2>(mesh,triangle)
              .accepted_conservatively())
        return false;
    }
    for(const GeometryVariable& variable : variables)
    {
      if(variable.kind != VariableKind::CurveParameter) continue;
      if(!planarP2BoundaryIsIntersectionFreeAroundPoint(
             mesh,variable.point))
        return false;
    }
    return true;
  }

  double globalObjective() const
  {
    bool invalid = false;
    double minimum = 0.;
    double maximum = 0.;
    double average = 0.;
    double value = 0.;
    if(objective == QuaFun::SizeShape)
      value = getmetquamesh<MetricFieldAnalytical,QuaFun::SizeShape>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
    else if(objective == QuaFun::StepDistance)
      value = getmetquamesh<MetricFieldAnalytical,QuaFun::StepDistance>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
    else
      throw std::runtime_error("Unsupported global smoothing objective");
    if(invalid || !std::isfinite(value))
      throw std::runtime_error("Invalid global objective value");
    return value;
  }

  bool evaluateValidValue(const std::vector<double>& x,
                          double& value,
                          bool updateBest)
  {
    try
    {
      apply(x);
      if(!geometryIsValid())
      {
        invalidValueEvaluations++;
        return false;
      }
      value = globalObjective();
      validValueEvaluations++;
      if(updateBest && value < bestValue)
      {
        bestValue = value;
        best = x;
      }
      return true;
    }
    catch(const std::exception&)
    {
      invalidValueEvaluations++;
      return false;
    }
  }

  double invalidPenalty(const std::vector<double>& x,
                        std::vector<double>* gradient) const
  {
    constexpr double base = 1.e12;
    constexpr double scale = 1.e6;
    double squaredDistance = 0.;
    for(std::size_t index = 0; index < x.size(); index++)
    {
      const double difference = x[index] - initial[index];
      squaredDistance += difference*difference;
      if(gradient != nullptr) (*gradient)[index] = 2.*scale*difference;
    }
    return base + scale*squaredDistance;
  }

  void frozenGradient(std::vector<double>& gradient)
  {
    std::vector<std::array<double,2>> physicalGradient(
        mesh.npoin,std::array<double,2>{0.,0.});
    const auto derivative = get_d_quafun<MetricFieldAnalytical,2,2>(objective);
    const int localNodeCount = getnnode(2,mesh.curdeg);
    for(int triangle = 0; triangle < mesh.nface; triangle++)
    {
      if(isdeadent(triangle,mesh.fac2poi)) continue;
      for(int localNode = 0; localNode < localNodeCount; localNode++)
      {
        double localGradient[2] = {0.,0.};
        (void)derivative(
            mesh,AsDeg::Pk,AsDeg::Pk,triangle,localNode,
            mesh.getBasis(),DifVar::None,localGradient,nullptr,1.0);
        const int point = mesh.fac2poi(triangle,localNode);
        physicalGradient[point][0] += localGradient[0];
        physicalGradient[point][1] += localGradient[1];
      }
    }

    for(std::size_t index = 0; index < variables.size(); index++)
    {
      const GeometryVariable& variable = variables[index];
      if(variable.kind == VariableKind::PhysicalX)
        gradient[index] = physicalGradient[variable.point][0];
      else if(variable.kind == VariableKind::PhysicalY)
        gradient[index] = physicalGradient[variable.point][1];
      else
      {
        double evaluation[18];
        const double parameter[2]
            = {mesh.bpo2rbi(variable.boundaryRecord,0),0.};
        const int status = EG_evaluate(variable.cadEdge,parameter,evaluation);
        if(status != EGADS_SUCCESS)
          throw std::runtime_error("EG_evaluate failed for CAD tangent");
        gradient[index]
            = physicalGradient[variable.point][0]*evaluation[3]
            + physicalGradient[variable.point][1]*evaluation[4];
      }
    }
  }

  void finiteDifferenceGradient(const std::vector<double>& x,
                                double baseValue,
                                std::vector<double>& gradient)
  {
    std::vector<double> perturbed = x;
    for(std::size_t index = 0; index < x.size(); index++)
    {
      const double interval = upper[index] - lower[index];
      double step = std::cbrt(std::numeric_limits<double>::epsilon())
                  * std::max({1.,std::abs(x[index]),interval});
      step = std::min(step,0.25*interval);
      const bool canGoLower = x[index] - step >= lower[index];
      const bool canGoUpper = x[index] + step <= upper[index];
      double lowerValue = 0.;
      double upperValue = 0.;
      bool lowerValid = false;
      bool upperValid = false;
      if(canGoLower)
      {
        perturbed[index] = x[index] - step;
        lowerValid = evaluateValidValue(perturbed,lowerValue,false);
      }
      if(canGoUpper)
      {
        perturbed[index] = x[index] + step;
        upperValid = evaluateValidValue(perturbed,upperValue,false);
      }
      perturbed[index] = x[index];

      if(lowerValid && upperValid)
        gradient[index] = (upperValue - lowerValue)/(2.*step);
      else if(upperValid)
        gradient[index] = (upperValue - baseValue)/step;
      else if(lowerValid)
        gradient[index] = (baseValue - lowerValue)/step;
      else
        gradient[index] = 0.;
    }
    apply(x);
  }

  double evaluate(const std::vector<double>& x,
                  std::vector<double>& gradient)
  {
    callbackCount++;
    double value = 0.;
    const bool valid = evaluateValidValue(x,value,true);
    if(!valid)
    {
      value = invalidPenalty(x,gradient.empty() ? nullptr : &gradient);
    }
    else if(!gradient.empty())
    {
      if(gradientMode == GradientMode::FiniteDifference)
        finiteDifferenceGradient(x,value,gradient);
      else
        frozenGradient(gradient);
    }

    double gradientNorm = 0.;
    for(const double entry : gradient) gradientNorm += entry*entry;
    gradientNorm = std::sqrt(gradientNorm);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - startTime).count();
    history << callbackCount << ',' << valid << ',' << value << ','
            << gradientNorm << ',' << validValueEvaluations << ','
            << invalidValueEvaluations << ',' << seconds << '\n';
    history.flush();
    std::cout << std::setprecision(17)
              << "Q2_GLOBAL_SMOOTHING {\"evaluation\":" << callbackCount
              << ",\"valid\":" << (valid ? "true" : "false")
              << ",\"objective_energy\":" << value
              << ",\"gradient_norm\":" << gradientNorm << "}\n";
    return value;
  }
};

double nloptObjective(const std::vector<double>& x,
                      std::vector<double>& gradient,
                      void* data)
{
  return static_cast<GlobalSmoothingProblem*>(data)->evaluate(x,gradient);
}

void constructVariables(GlobalSmoothingProblem& problem)
{
  Mesh<MetricFieldAnalytical>& mesh = problem.mesh;
  std::map<int,std::vector<std::size_t>> curveVariables;
  for(int point = 0; point < mesh.npoin; point++)
  {
    if(mesh.isdeadpoint(point)) continue;
    const int boundaryRecord = mesh.poi2bpo[point];
    if(boundaryRecord < 0)
    {
      problem.variables.push_back(
          {VariableKind::PhysicalX,point,-1,-1,nullptr});
      problem.initial.push_back(mesh.coord(point,0));
      problem.variables.push_back(
          {VariableKind::PhysicalY,point,-1,-1,nullptr});
      problem.initial.push_back(mesh.coord(point,1));
      continue;
    }

    const int cadDimension = mesh.bpo2ibi(boundaryRecord,1);
    if(cadDimension == 0) continue;
    if(cadDimension != 1)
      throw std::runtime_error(
          "The focused 2D test supports interior and CAD-curve points only");
    const int meshEdge = mesh.bpo2ibi(boundaryRecord,2);
    if(meshEdge < 0 || meshEdge >= mesh.nedge)
      throw std::runtime_error("Invalid owner edge for boundary point");
    const int cadReference = mesh.edg2ref[meshEdge];
    if(cadReference < 0 || cadReference >= mesh.CAD.ncaded)
      throw std::runtime_error("Invalid CAD edge reference");
    const std::size_t index = problem.variables.size();
    problem.variables.push_back({
        VariableKind::CurveParameter,point,boundaryRecord,cadReference,
        mesh.CAD.cad2edg[cadReference]});
    problem.initial.push_back(mesh.bpo2rbi(boundaryRecord,0));
    curveVariables[cadReference].push_back(index);
  }

  double minimumCoordinate[2] = {
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  double maximumCoordinate[2] = {
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for(int point = 0; point < mesh.npoin; point++)
  {
    if(mesh.isdeadpoint(point)) continue;
    for(int component = 0; component < 2; component++)
    {
      minimumCoordinate[component]
          = std::min(minimumCoordinate[component],mesh.coord(point,component));
      maximumCoordinate[component]
          = std::max(maximumCoordinate[component],mesh.coord(point,component));
    }
  }
  problem.lower.resize(problem.variables.size());
  problem.upper.resize(problem.variables.size());
  for(std::size_t index = 0; index < problem.variables.size(); index++)
  {
    const GeometryVariable& variable = problem.variables[index];
    if(variable.kind == VariableKind::PhysicalX)
    {
      problem.lower[index] = minimumCoordinate[0];
      problem.upper[index] = maximumCoordinate[0];
    }
    else if(variable.kind == VariableKind::PhysicalY)
    {
      problem.lower[index] = minimumCoordinate[1];
      problem.upper[index] = maximumCoordinate[1];
    }
  }

  for(auto& [cadReference,indices] : curveVariables)
  {
    ego cadEdge = mesh.CAD.cad2edg[cadReference];
    double range[4];
    int periodic = 0;
    const int status = EG_getRange(cadEdge,range,&periodic);
    if(status != EGADS_SUCCESS || periodic != 0)
      throw std::runtime_error(
          "The focused test requires non-periodic CAD boundary curves");
    std::sort(indices.begin(),indices.end(),[&](std::size_t first,
                                               std::size_t second)
    {
      return problem.initial[first] < problem.initial[second];
    });
    const double edgeLength = range[1] - range[0];
    const double margin = 256.*std::numeric_limits<double>::epsilon()
                        * std::max(1.,std::abs(edgeLength));
    for(std::size_t position = 0; position < indices.size(); position++)
    {
      const std::size_t index = indices[position];
      const double previous = position == 0
                            ? range[0]
                            : problem.initial[indices[position - 1]];
      const double next = position + 1 == indices.size()
                        ? range[1]
                        : problem.initial[indices[position + 1]];
      problem.lower[index] = 0.5*(previous + problem.initial[index]) + margin;
      problem.upper[index] = 0.5*(problem.initial[index] + next) - margin;
      if(!(problem.lower[index] < problem.initial[index]
           && problem.initial[index] < problem.upper[index]))
        throw std::runtime_error("Invalid CAD parameter trust interval");
    }
  }
  problem.best = problem.initial;
}

std::pair<double,double> directionalDerivativeAudit(
    GlobalSmoothingProblem& problem)
{
  std::vector<double> direction(problem.initial.size());
  for(std::size_t index = 0; index < direction.size(); index++)
    direction[index] = std::sin(static_cast<double>(index + 1));
  const double norm = std::sqrt(std::inner_product(
      direction.begin(),direction.end(),direction.begin(),0.));
  for(double& entry : direction) entry /= norm;

  problem.apply(problem.initial);
  std::vector<double> frozen(problem.initial.size(),0.);
  problem.frozenGradient(frozen);
  const double frozenDirectional = std::inner_product(
      frozen.begin(),frozen.end(),direction.begin(),0.);

  double maximumStep = std::numeric_limits<double>::infinity();
  for(std::size_t index = 0; index < direction.size(); index++)
  {
    if(direction[index] > 0.)
      maximumStep = std::min(
          maximumStep,
          (problem.upper[index] - problem.initial[index])/direction[index]);
    else if(direction[index] < 0.)
      maximumStep = std::min(
          maximumStep,
          (problem.lower[index] - problem.initial[index])/direction[index]);
  }
  const double step = std::min(1.e-6,0.1*maximumStep);
  std::vector<double> lowerPoint = problem.initial;
  std::vector<double> upperPoint = problem.initial;
  for(std::size_t index = 0; index < direction.size(); index++)
  {
    lowerPoint[index] -= step*direction[index];
    upperPoint[index] += step*direction[index];
  }
  double lowerValue = 0.;
  double upperValue = 0.;
  if(!problem.evaluateValidValue(lowerPoint,lowerValue,false)
     || !problem.evaluateValidValue(upperPoint,upperValue,false))
    throw std::runtime_error("Directional gradient audit left valid geometry");
  problem.apply(problem.initial);
  return {frozenDirectional,(upperValue - lowerValue)/(2.*step)};
}

} // namespace

int main(int argc, char** argv)
{
  try
  {
    if(argc == 1)
    {
      std::cerr
          << "Usage: smoothQ2GlobalNodes --in input.meshb --cad model.egads "
             "--anamet ID --sclmet SCALE --tardeg 2 --out final.meshb "
             "--elevated-output elevated.meshb --history-output history.csv "
             "[--smoothing-objective sizeshape|step-distance-sv] "
             "[--gradient-mode finite-difference|frozen] "
             "[--maximum-evaluations 100] [Metris options]\n";
      return 2;
    }

    CustomOptions options = extractCustomOptions(argc,argv);
    for(const auto& output : {options.elevatedOutput,options.historyOutput})
    {
      if(std::filesystem::exists(output))
        throw std::runtime_error("Refusing to overwrite: " + output.string());
      if(!output.parent_path().empty())
        std::filesystem::create_directories(output.parent_path());
    }

    MetrisRunner run(
        static_cast<int>(options.metrisArgv.size()),
        options.metrisArgv.data());
    if(run.metricFE)
      throw std::runtime_error("This focused test requires an analytic metric");

    auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*run.msh_g);
    if(mesh.idim != 2 || mesh.get_tdim() != 2 || mesh.curdeg != 1)
      throw std::runtime_error("Input must be a full-dimensional planar P1 mesh");
    if(run.param->usrTarDeg != 2)
      throw std::runtime_error("This focused test requires --tardeg 2");
    if(options.objective == QuaFun::StepDistance
       && !run.param->step_distance_shape_volume)
      throw std::runtime_error(
          "step-distance-sv requires --step-distance-shape-volume");

    const int originalPointCount = mesh.npoin;
    run.degElevate();
    mesh.setBasis(FEBasis::Lagrange);
    if(mesh.curdeg != 2)
      throw std::runtime_error("Degree elevation did not produce P2 geometry");
    writeMesh(options.elevatedOutput.string(),mesh);

    std::ofstream history(options.historyOutput);
    if(!history.good())
      throw std::runtime_error(
          "Cannot open history: " + options.historyOutput.string());
    history << std::setprecision(17)
            << "callback,valid,objective_energy,gradient_norm,"
               "valid_value_evaluations,invalid_value_evaluations,seconds\n";

    GlobalSmoothingProblem problem{
        mesh,options.objective,options.gradientMode,{}, {}, {}, {}, {},
        std::numeric_limits<double>::infinity(),history,
        std::chrono::steady_clock::now()};
    constructVariables(problem);
    double initialValue = 0.;
    if(!problem.evaluateValidValue(problem.initial,initialValue,true))
      throw std::runtime_error("Initial elevated mesh is invalid");

    const auto [frozenDirectional,finiteDifferenceDirectional]
        = directionalDerivativeAudit(problem);
    const double directionalDifference
        = std::abs(frozenDirectional - finiteDifferenceDirectional);
    const double directionalScale = std::max(
        {1.,std::abs(frozenDirectional),
         std::abs(finiteDifferenceDirectional)});

    std::cout << std::setprecision(17)
              << "Q2_GLOBAL_SETUP {\"objective\":\""
              << options.objectiveName << "\",\"gradient\":\""
              << options.gradientName << "\",\"p1_points\":"
              << originalPointCount << ",\"q2_points\":" << mesh.npoin
              << ",\"variables\":" << problem.variables.size()
              << ",\"initial_energy\":" << initialValue
              << ",\"frozen_directional_derivative\":"
              << frozenDirectional
              << ",\"finite_difference_directional_derivative\":"
              << finiteDifferenceDirectional
              << ",\"directional_relative_mismatch\":"
              << directionalDifference/directionalScale << "}\n";

    nlopt::opt optimizer(nlopt::LD_LBFGS,problem.variables.size());
    optimizer.set_lower_bounds(problem.lower);
    optimizer.set_upper_bounds(problem.upper);
    optimizer.set_min_objective(nloptObjective,&problem);
    optimizer.set_ftol_rel(1.e-10);
    optimizer.set_xtol_rel(1.e-8);
    optimizer.set_maxeval(options.maximumEvaluations);

    std::vector<double> optimized = problem.initial;
    double optimizerValue = initialValue;
    nlopt::result result = nlopt::FAILURE;
    try
    {
      result = optimizer.optimize(optimized,optimizerValue);
    }
    catch(const nlopt::roundoff_limited&)
    {
      std::cerr << "NLopt stopped because further progress was roundoff-limited\n";
    }
    catch(const std::exception& exception)
    {
      // Some NLopt algorithms report a failed line search as an exception.
      // The callback has independently retained the best valid state, so this
      // is a recoverable optimizer termination rather than a failed test.
      std::cerr << "NLopt stopped without a success code: "
                << exception.what() << '\n';
    }

    problem.apply(problem.best);
    if(!problem.geometryIsValid())
      throw std::runtime_error("Best global-smoothing state is invalid");
    const double finalValue = problem.globalObjective();
    if(finalValue > initialValue
       + 1.e-12*std::max(1.,std::abs(initialValue)))
      throw std::runtime_error("Global smoothing increased the objective");

    std::cout << std::setprecision(17)
              << "Q2_GLOBAL_RESULT {\"nlopt_result\":"
              << static_cast<int>(result)
              << ",\"callbacks\":" << problem.callbackCount
              << ",\"valid_value_evaluations\":"
              << problem.validValueEvaluations
              << ",\"invalid_value_evaluations\":"
              << problem.invalidValueEvaluations
              << ",\"initial_energy\":" << initialValue
              << ",\"final_energy\":" << finalValue
              << ",\"relative_reduction\":"
              << (initialValue - finalValue)/initialValue << "}\n";

    run.writeOutputs();
  }
  catch(const std::exception& exception)
  {
    std::cerr << "smoothQ2GlobalNodes: " << exception.what() << '\n';
    return 1;
  }
  return 0;
}

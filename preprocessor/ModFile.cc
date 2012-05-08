/*
 * Copyright (C) 2006-2012 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <typeinfo>
#ifndef _WIN32
# include <unistd.h>
#endif

#include "ModFile.hh"
#include "ConfigFile.hh"

ModFile::ModFile() : expressions_tree(symbol_table, num_constants, external_functions_table),
                     dynamic_model(symbol_table, num_constants, external_functions_table),
                     trend_dynamic_model(symbol_table, num_constants, external_functions_table),
                     static_model(symbol_table, num_constants, external_functions_table),
                     steady_state_model(symbol_table, num_constants, external_functions_table, static_model),
                     linear(false), block(false), byte_code(false),
                     use_dll(false), no_static(false), nonstationary_variables(false)
{
}

ModFile::~ModFile()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    delete (*it);
}

void
ModFile::evalAllExpressions(bool warn_uninit)
{
  cout << "Evaluating expressions...";

  // Loop over all statements, and fill global eval context if relevant
  for (vector<Statement *>::const_iterator it = statements.begin(); it != statements.end(); it++)
    {
      InitParamStatement *ips = dynamic_cast<InitParamStatement *>(*it);
      if (ips)
        ips->fillEvalContext(global_eval_context);

      InitOrEndValStatement *ies = dynamic_cast<InitOrEndValStatement *>(*it);
      if (ies)
        ies->fillEvalContext(global_eval_context);

      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass)
        lpass->fillEvalContext(global_eval_context);
    }

  // Evaluate model local variables
  dynamic_model.fillEvalContext(global_eval_context);

  cout << "done" << endl;

  // Check if some symbols are not initialized, and give them a zero value then
  for (int id = 0; id <= symbol_table.maxID(); id++)
    {
      SymbolType type = symbol_table.getType(id);
      if ((type == eEndogenous || type == eExogenous || type == eExogenousDet
           || type == eParameter || type == eModelLocalVariable)
          && global_eval_context.find(id) == global_eval_context.end())
        {
          if (warn_uninit)
            cerr << "WARNING: can't find a numeric initial value for " << symbol_table.getName(id) << ", using zero" << endl;
          global_eval_context[id] = 0;
        }
    }
}

void
ModFile::addStatement(Statement *st)
{
  statements.push_back(st);
}

void
ModFile::addStatementAtFront(Statement *st)
{
  statements.insert(statements.begin(), st);
}

void
ModFile::checkPass()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->checkPass(mod_file_struct);

  // Check the steady state block
  steady_state_model.checkPass(mod_file_struct.ramsey_policy_present);

  // If order option has not been set, default to 2
  if (!mod_file_struct.order_option)
    mod_file_struct.order_option = 2;

  bool stochastic_statement_present = mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present
    || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present;

  // Allow empty model only when doing a standalone BVAR estimation
  if (dynamic_model.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.simul_present
          || stochastic_statement_present))
    {
      cerr << "ERROR: At least one model equation must be declared!" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.simul_present && stochastic_statement_present)
    {
      cerr << "ERROR: A .mod file cannot contain both a simul command and one of {stoch_simul, estimation, osr, ramsey_policy}" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.k_order_solver && byte_code)
    {
      cerr << "ERROR: 'k_order_solver' (which is implicit if order >= 3), is not yet compatible with 'bytecode'." << endl;
      exit(EXIT_FAILURE);
    }

  if (use_dll && (block || byte_code))
    {
      cerr << "ERROR: In 'model' block, 'use_dll' option is not compatible with 'block' or 'bytecode'" << endl;
      exit(EXIT_FAILURE);
    }
  if (block || byte_code)
    if (dynamic_model.isModelLocalVariableUsed())
      {
        cerr << "ERROR: In 'model' block, 'block' or 'bytecode' options are not yet compatible with pound expressions" << endl;
        exit(EXIT_FAILURE);
      }
  if ((stochastic_statement_present || mod_file_struct.check_present || mod_file_struct.steady_present) && no_static)
    {
      cerr << "ERROR: no_static option is incompatible with stoch_simul, estimation, osr, ramsey_policy, steady and check commands" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.dsge_var_estimated)
    if (!mod_file_struct.dsge_prior_weight_in_estimated_params)
      {
        cerr << "ERROR: When estimating a DSGE-VAR model and estimating the weight of the prior, dsge_prior_weight must "
             << "be referenced in the estimated_params block." << endl;
        exit(EXIT_FAILURE);
      }

  if (symbol_table.exists("dsge_prior_weight"))
    {
      if (symbol_table.getType("dsge_prior_weight") != eParameter)
        {
          cerr << "ERROR: dsge_prior_weight may only be used as a parameter." << endl;
          exit(EXIT_FAILURE);
        }
      else
        cerr << "WARNING: When estimating a DSGE-Var, declaring dsge_prior_weight as a parameter is deprecated. "
             <<  "The preferred method is to do this via the dsge_var option in the estimation statement." << endl;

      if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
        {
          cerr << "ERROR: dsge_prior_weight can either be declared as a parameter (deprecated) or via the dsge_var option "
               << "to the estimation statement (preferred), but not both." << endl;
          exit(EXIT_FAILURE);
        }

      if (!mod_file_struct.dsge_prior_weight_initialized && !mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: If dsge_prior_weight is declared as a parameter, it must either be initialized or placed in the "
               << "estimated_params block." << endl;
          exit(EXIT_FAILURE);
        }

      if (mod_file_struct.dsge_prior_weight_initialized && mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: dsge_prior_weight cannot be both initalized and estimated." << endl;
          exit(EXIT_FAILURE);
        }
    }

  if (mod_file_struct.dsge_prior_weight_in_estimated_params)
    if (!mod_file_struct.dsge_var_estimated && !mod_file_struct.dsge_var_calibrated.empty())
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, the prior weight cannot be calibrated "
             << "via the dsge_var option in the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }
    else if (!mod_file_struct.dsge_var_estimated && !symbol_table.exists("dsge_prior_weight"))
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, it must either be declared as a parameter "
             << "(deprecated) or the dsge_var option must be passed to the estimation statement (preferred)." << endl;
        exit(EXIT_FAILURE);
      }
}

void
ModFile::transformPass()
{
  if (symbol_table.predeterminedNbr() > 0)
    dynamic_model.transformPredeterminedVariables();

  // Create auxiliary vars for Expectation operator
  dynamic_model.substituteExpectation(mod_file_struct.partial_information);

  if (nonstationary_variables)
    {
      dynamic_model.detrendEquations();
      dynamic_model.cloneDynamic(trend_dynamic_model);
      dynamic_model.removeTrendVariableFromEquations();
    }

  if (mod_file_struct.stoch_simul_present
      || mod_file_struct.estimation_present
      || mod_file_struct.osr_present
      || mod_file_struct.ramsey_policy_present)
    {
      // In stochastic models, create auxiliary vars for leads and lags greater than 2, on both endos and exos
      dynamic_model.substituteEndoLeadGreaterThanTwo(false);
      dynamic_model.substituteExoLead(false);
      dynamic_model.substituteEndoLagGreaterThanTwo(false);
      dynamic_model.substituteExoLag(false);
    }
  else
    {
      // In deterministic models, create auxiliary vars for leads and lags endogenous greater than 2, only on endos (useless on exos)
      dynamic_model.substituteEndoLeadGreaterThanTwo(true);
      dynamic_model.substituteEndoLagGreaterThanTwo(true);
    }

  if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
    try
      {
        int sid = symbol_table.addSymbol("dsge_prior_weight", eParameter);
        if (!mod_file_struct.dsge_var_calibrated.empty())
          addStatementAtFront(new InitParamStatement(sid,
                                                     expressions_tree.AddNonNegativeConstant(mod_file_struct.dsge_var_calibrated),
                                                     symbol_table));
      }
    catch (SymbolTable::AlreadyDeclaredException &e)
      {
        cerr << "ERROR: dsge_prior_weight should not be declared as a model variable / parameter "
             << "when the dsge_var option is passed to the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }

  // Freeze the symbol table
  symbol_table.freeze();

  /*
    Enforce the same number of equations and endogenous, except in two cases:
    - ramsey_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
  */
  if (!mod_file_struct.ramsey_policy_present
      && !(mod_file_struct.bvar_present && dynamic_model.equation_number() == 0)
      && (dynamic_model.equation_number() != symbol_table.endo_nbr()))
    {
      cerr << "ERROR: There are " << dynamic_model.equation_number() << " equations but " << symbol_table.endo_nbr() << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }
  
  if (symbol_table.exo_det_nbr() > 0 && mod_file_struct.simul_present)
    {
      cerr << "ERROR: A .mod file cannot contain both a simul command and varexo_det declaration (all exogenous variables are deterministic in this case)" << endl;
      exit(EXIT_FAILURE);
    }
  
  cout << "Found " << dynamic_model.equation_number() << " equation(s)." << endl;

  if (symbol_table.exists("dsge_prior_weight"))
    if (mod_file_struct.bayesian_irf_present)
      {
        if (symbol_table.exo_nbr() != symbol_table.observedVariablesNbr())
          {
            cerr << "ERROR: When estimating a DSGE-Var and the bayesian_irf option is passed to the estimation "
                 << "statement, the number of shocks must equal the number of observed variables." << endl;
            exit(EXIT_FAILURE);
          }
      }
    else
      if (symbol_table.exo_nbr() < symbol_table.observedVariablesNbr())
        {
          cerr << "ERROR: When estimating a DSGE-Var, the number of shocks must be "
               << "greater than or equal to the number of observed variables." << endl;
          exit(EXIT_FAILURE);
        }
}

void
ModFile::computingPass(bool no_tmp_terms)
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  bool dynamic_model_needed = mod_file_struct.simul_present || mod_file_struct.check_present || mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present || mod_file_struct.identification_present;
  if (dynamic_model.equation_number() > 0)
    {
      if (nonstationary_variables)
        trend_dynamic_model.runTrendTest(global_eval_context);

      // Compute static model and its derivatives
      dynamic_model.toStatic(static_model);
      if (!no_static)
        {
          static_model.initializeVariablesAndEquations();
          static_model.computingPass(global_eval_context, no_tmp_terms, false, block, byte_code);
        }
      // Set things to compute for dynamic model
      if (dynamic_model_needed)
        {
          dynamic_model.initializeVariablesAndEquations();
          if (mod_file_struct.simul_present)
            dynamic_model.computingPass(true, false, false, false, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
          else
            {
              if (mod_file_struct.order_option < 1 || mod_file_struct.order_option > 3)
                {
                  cerr << "ERROR: Incorrect order option..." << endl;
                  exit(EXIT_FAILURE);
                }
              bool hessian = mod_file_struct.order_option >= 2 || mod_file_struct.identification_present;
              bool thirdDerivatives = mod_file_struct.order_option == 3;
              bool paramsDerivatives = mod_file_struct.identification_present;
              dynamic_model.computingPass(true, hessian, thirdDerivatives, paramsDerivatives, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
            }
        }
      else
        dynamic_model.computingPass(true, true, false, false, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
    }

  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->computingPass();
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all, bool console, const ConfigFile &config_file
#if defined(_WIN32) || defined(__CYGWIN32__)
                          , bool cygwin, bool msvc
#endif
                          ) const
{
  ofstream mOutputFile;
  bool dynamic_model_needed = mod_file_struct.simul_present || mod_file_struct.check_present || mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present || mod_file_struct.identification_present;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : main Dynare file " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (clear_all)
    {
      mOutputFile << "clear all" << endl
	// this is a work-around for a bug in Octave 3.2
		  << "clear global" << endl;
    }

  mOutputFile << "tic;" << endl
              << "global M_ oo_ options_ ys0_ ex0_" << endl
              << "options_ = [];" << endl
              << "M_.fname = '" << basename << "';" << endl
              << "%" << endl
              << "% Some global variables initialization" << endl
              << "%" << endl
              << "global_initialization;" << endl
              << "diary off;" << endl
              << "logname_ = '" << basename << ".log';" << endl
              << "if exist(logname_, 'file')" << endl
              << "    delete(logname_)" << endl
              << "end" << endl
              << "diary(logname_)" << endl;

  if (console)
    mOutputFile << "options_.console_mode = 1;" << endl;

  cout << "Processing outputs ...";

  symbol_table.writeOutput(mOutputFile);

  // Initialize M_.Sigma_e and M_.H
  mOutputFile << "M_.Sigma_e = zeros(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl;

  if (mod_file_struct.calibrated_measurement_errors)
    mOutputFile << "M_.H = zeros(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl;
  else
    mOutputFile << "M_.H = 0;" << endl;

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  mOutputFile << "options_.block=" << block << ";" << endl
              << "options_.bytecode=" << byte_code << ";" << endl
              << "options_.use_dll=" << use_dll << ";" << endl;

  config_file.writeCluster(mOutputFile);

  if (byte_code)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  // Erase possible remnants of previous runs
  if (block || byte_code || use_dll)
    mOutputFile << "if exist('" << basename << "_dynamic.m', 'file')" << endl
                << "    delete('" << basename << "_dynamic.m');" << endl
                << "end" << endl;

  if (byte_code)
    mOutputFile << "if exist('" << basename << "_static.m', 'file')" << endl
                << "    delete('" << basename << "_static.m');" << endl
                << "end" << endl;

  if (!use_dll)
    mOutputFile << "erase_compiled_function('" + basename + "_dynamic');" << endl;

  // Erase generated steady state file (see ticket #224)
  string steadystatefile = basename + "_steadystate.m";
  ifstream in(steadystatefile.c_str(), ios::binary);
  if (!in.fail())
    {
      string line;
      getline(in, line);
      if (!line.compare(STEADY_STATE_GENERATED_HEADER))
        {
          in.close();
          unlink(steadystatefile.c_str());
        }
      else
        in.close();
    }
  
#if defined(_WIN32) || defined(__CYGWIN32__)
  // If using USE_DLL with MSVC, check that the user didn't use a function not supported by MSVC (because MSVC doesn't comply with C99 standard)
  if (use_dll && msvc)
    {
      if (dynamic_model.isUnaryOpUsed(oAcosh))
        {
          cerr << "ERROR: acosh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAsinh))
        {
          cerr << "ERROR: asinh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAtanh))
        {
          cerr << "ERROR: atanh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isTrinaryOpUsed(oNormcdf))
        {
          cerr << "ERROR: normcdf() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
    }
#endif

  // Compile the dynamic MEX file for use_dll option
  if (use_dll)
    {
      mOutputFile << "if ~exist('OCTAVE_VERSION')" << endl;
      // Some mex commands are enclosed in an eval(), because otherwise it will make Octave fail
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (msvc)
        mOutputFile << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Dynamic\" " << basename << "_dynamic.c')" << endl;                                                                                                                                                                                                                                                                                                                                                                                   // MATLAB/Windows + Microsoft Visual C++
      else if (cygwin)
        mOutputFile << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Dynamic >> mex.def\" " << basename << "_dynamic.c')" << endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            // MATLAB/Windows + Cygwin g++
      else
        mOutputFile << "    error('When using the USE_DLL option, you must give either ''cygwin'' or ''msvc'' option to the ''dynare'' command')" << endl;
#else
# ifdef __linux__
      mOutputFile << "    eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_dynamic.c')" << endl; // MATLAB/Linux
# else // MacOS
      mOutputFile << "    eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' " << basename << "_dynamic.c')" << endl; // MATLAB/MacOS
# endif
#endif
      mOutputFile << "else" << endl // Octave
                  << "    if ~octave_ver_less_than('3.2.0')" << endl // Workaround for bug in Octave >= 3.2, see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=550823
                  << "        sleep(2)" << endl
                  << "    end" << endl
                  << "    mex " << basename << "_dynamic.c" << endl
                  << "end" << endl;
    }

  // Add path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "addpath " << basename << ";" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeOutput(mOutputFile, basename, block, byte_code, use_dll, mod_file_struct.order_option);
      if (!no_static)
        static_model.writeOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
    {
      (*it)->writeOutput(mOutputFile, basename);

      // Special treatment for initval block: insert initial values for the auxiliary variables
      InitValStatement *ivs = dynamic_cast<InitValStatement *>(*it);
      if (ivs != NULL)
        {
          static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
          ivs->writeOutputPostInit(mOutputFile);
        }

      // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass && !no_static)
        static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
    }

  // Remove path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "rmpath " << basename << ";" << endl;

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl;

  config_file.writeEndParallel(mOutputFile);

  mOutputFile << "diary off" << endl
              << endl << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;

  mOutputFile.close();

  // Create static and dynamic files
  if (dynamic_model.equation_number() > 0)
    {
      if (!no_static)
        static_model.writeStaticFile(basename, block, byte_code);

        dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);
        dynamic_model.writeParamsDerivativesFile(basename);
    }

  // Create steady state file
  steady_state_model.writeSteadyStateFile(basename, mod_file_struct.ramsey_policy_present);

  cout << "done" << endl;
}

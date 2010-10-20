/*
 * Copyright (C) 2003-2009 Dynare Team
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

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "NumericalInitialization.hh"

InitParamStatement::InitParamStatement(int symb_id_arg,
                                       const NodeID param_value_arg,
                                       const SymbolTable &symbol_table_arg) :
  symb_id(symb_id_arg),
  param_value(param_value_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitParamStatement::writeOutput(ostream &output, const string &basename) const
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << "M_.params( " << id << " ) = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  output << symbol_table.getName(symb_id) << " = M_.params( " << id << " );\n";
}

void
InitParamStatement::fillEvalContext(eval_context_type &eval_context) const
{
  try
    {
      eval_context[symb_id] = param_value->eval(eval_context);
    }
  catch (ExprNode::EvalException &e)
    {
      // Do nothing
    }
}

InitOrEndValStatement::InitOrEndValStatement(const init_values_type &init_values_arg,
                                             const SymbolTable &symbol_table_arg) :
  init_values(init_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitOrEndValStatement::fillEvalContext(eval_context_type &eval_context) const
{
  for (init_values_type::const_iterator it = init_values.begin();
       it != init_values.end(); it++)
    {
      try
        {
          eval_context[it->first] = (it->second)->eval(eval_context);
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for (init_values_type::const_iterator it = init_values.begin();
       it != init_values.end(); it++)
    {
      const int symb_id = it->first;
      const NodeID expression = it->second;

      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      if (type == eEndogenous)
        output << "oo_.steady_state";
      else if (type == eExogenous)
        output << "oo_.exo_steady_state";
      else if (type == eExogenousDet)
        output << "oo_.exo_det_steady_state";

      output << "( " << tsid << " ) = ";
      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitValStatement::InitValStatement(const init_values_type &init_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}

void
InitValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% INITVAL instructions" << endl
         << "%" << endl;
  // Writing initval block to set initial values for variables
  output << "options_.initval_file = 0;" << endl;

  writeInitValues(output);
}

void
InitValStatement::writeOutputPostInit(ostream &output) const
{
  output << "oo_.endo_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];" << endl
         << "if M_.exo_nbr > 0;" << endl
         << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];" << endl
         <<"end;" << endl
         << "if M_.exo_det_nbr > 0;" << endl
         << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];" << endl
         <<"end;" << endl;
}

EndValStatement::EndValStatement(const init_values_type &init_values_arg,
                                 const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}

void
EndValStatement::checkPass(ModFileStructure &mod_file_struct)
{
  if (mod_file_struct.shocks_present)
    {
      cerr << "ERROR: Putting a \"shocks\" block before an \"endval\" block is not permitted. Please swap the two blocks. This limitation will be removed in the next major release of Dynare." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EndValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% ENDVAL instructions" << endl
         << "%" << endl;
  // Writing endval block to set terminal values for variables
  output << "ys0_= oo_.steady_state;" << endl
         << "ex0_ = oo_.exo_steady_state;" << endl;

  writeInitValues(output);
}

HistValStatement::HistValStatement(const hist_values_type &hist_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  hist_values(hist_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HistValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% HISTVAL instructions" << endl
         << "%" << endl
         << "oo_.endo_simul = zeros(M_.endo_nbr,M_.maximum_lag);" << endl;

  for (hist_values_type::const_iterator it = hist_values.begin();
       it != hist_values.end(); it++)
    {
      int symb_id = it->first.first;
      int lag = it->first.second;
      const NodeID expression = it->second;

      SymbolType type = symbol_table.getType(symb_id);
      if ((type == eEndogenous && lag < 0) || type == eExogenous)
	{
	  const int new_symb_id = symbol_table.searchAuxiliaryVars(symb_id,lag);
	  if (new_symb_id != -1)
	    {
	      symb_id = new_symb_id;
	      lag = 0;
        type = eEndogenous;
	    }
	  else if (symbol_table.AuxVarsSize() > 0)
	    {
	      cerr << "Histval: this variable doesn't exist with such a lag in the model" << endl;
	      exit(EXIT_FAILURE);
	    }
	    
	}
      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      if (type == eEndogenous)
        output << "oo_.endo_simul( " << tsid << ", M_.maximum_lag + " << lag << ") = ";
      else if (type == eExogenous)
        output << "oo_.exo_simul( M_.maximum_lag + " << lag << ", " << tsid << " ) = ";
      else if (type != eExogenousDet)
        output << "oo_.exo_det_simul( M_.maximum_lag + " << lag  << ", " << tsid << " ) = ";

      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitvalFileStatement::InitvalFileStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
InitvalFileStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% INITVAL_FILE statement" << endl
         << "%" << endl
         << "options_.initval_file = 1;" << endl
         << "initvalf('" << filename << "');" << endl;
}

HomotopyStatement::HomotopyStatement(const homotopy_values_type &homotopy_values_arg,
                                     const SymbolTable &symbol_table_arg) :
  homotopy_values(homotopy_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HomotopyStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% HOMOTOPY_SETUP instructions" << endl
         << "%" << endl
         << "options_.homotopy_values = [];" << endl;

  for (homotopy_values_type::const_iterator it = homotopy_values.begin();
       it != homotopy_values.end(); it++)
    {
      const int &symb_id = it->first;
      const NodeID expression1 = it->second.first;
      const NodeID expression2 = it->second.second;

      const SymbolType type = symbol_table.getType(symb_id);
      const int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      output << "options_.homotopy_values = vertcat(options_.homotopy_values, [ " << type << ", " << tsid << ", ";
      if (expression1 != NULL)
        expression1->writeOutput(output);
      else
        output << "NaN";
      output << ", ";
      expression2->writeOutput(output);
      output << "]);" << endl;
    }
}

SaveParamsAndSteadyStateStatement::SaveParamsAndSteadyStateStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
SaveParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "save_params_and_steady_state('" << filename << "');" << endl;
}

LoadParamsAndSteadyStateStatement::LoadParamsAndSteadyStateStatement(const string &filename,
                                                                     const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
  cout << "Reading " << filename << "." << endl;

  ifstream f;
  f.open(filename.c_str(), ios::in);
  if (f.fail())
    {
      cerr << "ERROR: Can't open " << filename << endl;
      exit(EXIT_FAILURE);
    }

  while (true)
    {
      string symb_name, value;
      f >> symb_name >> value;
      if (f.eof())
        break;

      try
        {
          int symb_id = symbol_table.getID(symb_name);
          content[symb_id] = value;
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          cerr << "WARNING: Unknown symbol " << symb_name << " in " << filename << endl;
        }
    }
}

void
LoadParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename) const
{
  for (map<int, string>::const_iterator it = content.begin();
       it != content.end(); it++)
    {
      switch (symbol_table.getType(it->first))
        {
        case eParameter:
          output << "M_.params";
          break;
        case eEndogenous:
          output << "oo_.steady_state";
          break;
        case eExogenous:
          output << "oo_.exo_steady_state";
          break;
        case eExogenousDet:
          output << "oo_.exo_det_steady_state";
          break;
        default:
          cerr << "ERROR: Unsupported variable type for " << symbol_table.getName(it->first) << " in load_params_and_steady_state" << endl;
          exit(EXIT_FAILURE);
        }

      int tsid = symbol_table.getTypeSpecificID(it->first) + 1;
      output << "(" << tsid << ") = " << it->second << ";" << endl;
    }
}

void
LoadParamsAndSteadyStateStatement::fillEvalContext(eval_context_type &eval_context) const
{
  for (map<int, string>::const_iterator it = content.begin();
       it != content.end(); it++)
    eval_context[it->first] = atof(it->second.c_str());
}

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
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cerrno>
#include <algorithm>
#include "DynamicModel.hh"

// For mkdir() and chdir()
#ifdef _WIN32
# include <direct.h>
#else
# include <unistd.h>
# include <sys/stat.h>
# include <sys/types.h>
#endif

DynamicModel::DynamicModel(SymbolTable &symbol_table_arg,
                           NumericalConstants &num_constants_arg) :
  ModelTree(symbol_table_arg, num_constants_arg),
  max_lag(0), max_lead(0),
  max_endo_lag(0), max_endo_lead(0),
  max_exo_lag(0), max_exo_lead(0),
  max_exo_det_lag(0), max_exo_det_lead(0),
  dynJacobianColsNbr(0),
  global_temporary_terms(true),
  cutoff(1e-15),
  mfs(0)
{
}

VariableNode *
DynamicModel::AddVariable(int symb_id, int lag)
{
  return AddVariableInternal(symb_id, lag);
}

void
DynamicModel::compileDerivative(ofstream &code_file, int eq, int symb_id, int lag, map_idx_type &map_idx) const
{
  first_derivatives_type::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, symb_id), lag)));
  if (it != first_derivatives.end())
    (it->second)->compile(code_file, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file);
    }
}

void
DynamicModel::compileChainRuleDerivative(ofstream &code_file, int eqr, int varr, int lag, map_idx_type &map_idx) const
{
  map<pair<int, pair<int, int> >, NodeID>::const_iterator it = first_chain_rule_derivatives.find(make_pair(eqr, make_pair(varr, lag)));
  if (it != first_chain_rule_derivatives.end())
    (it->second)->compile(code_file, false, temporary_terms, map_idx, true, false);
  else
    {
      FLDZ_ fldz;
      fldz.write(code_file);
    }
}

void
DynamicModel::computeTemporaryTermsOrdered()
{
  map<NodeID, pair<int, int> > first_occurence;
  map<NodeID, int> reference_count;
  BinaryOpNode *eq_node;
  first_derivatives_type::const_iterator it;
  first_chain_rule_derivatives_type::const_iterator it_chr;
  ostringstream tmp_s;
  v_temporary_terms.clear();
  map_idx.clear();

  unsigned int nb_blocks = getNbBlocks();
  v_temporary_terms = vector<vector<temporary_terms_type> >(nb_blocks);
  v_temporary_terms_inuse = vector<temporary_terms_inuse_type>(nb_blocks);
  temporary_terms.clear();

  if (!global_temporary_terms)
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          reference_count.clear();
          temporary_terms.clear();
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_type>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
            }
          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);
          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  block_size-1);

          set<int> temporary_terms_in_use;
          temporary_terms_in_use.clear();
          v_temporary_terms_inuse[block] = temporary_terms_in_use;
        }
    }
  else
    {
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Compute the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          v_temporary_terms[block] = vector<temporary_terms_type>(block_size);
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms,  i);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, i);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
            }
          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, block, v_temporary_terms, block_size-1);
        }
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          // Collect the temporary terms reordered
          unsigned int block_size = getBlockSize(block);
          unsigned int block_nb_mfs = getBlockMfs(block);
          unsigned int block_nb_recursives = block_size - block_nb_mfs;
          set<int> temporary_terms_in_use;
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i < block_nb_recursives && isBlockEquationRenormalized(block, i))
                getBlockEquationRenormalizedNodeID(block, i)->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
              else
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  eq_node->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
                }
            }
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              NodeID id = it->second.second;
              id->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
            }
          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            it->second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            it->second->collectTemporary_terms(temporary_terms, temporary_terms_in_use, block);
          v_temporary_terms_inuse[block] = temporary_terms_in_use;
        }
      // Add a mapping form node ID to temporary terms order
      int j = 0;
      for (temporary_terms_type::const_iterator it = temporary_terms.begin();
           it != temporary_terms.end(); it++)
        map_idx[(*it)->idx] = j++;
    }
}

void
DynamicModel::writeModelEquationsOrdered_M(const string &dynamic_basename) const
{
  string tmp_s, sps;
  ostringstream tmp_output, tmp1_output, global_output;
  NodeID lhs = NULL, rhs = NULL;
  BinaryOpNode *eq_node;
  ostringstream Uf[symbol_table.endo_nbr()];
  map<NodeID, int> reference_count;
  temporary_terms_type local_temporary_terms;
  ofstream  output;
  int nze, nze_exo, nze_other_endo;
  vector<int> feedback_variables;
  ExprNodeOutputType local_output_type;

  if (global_temporary_terms)
    {
      local_output_type = oMatlabDynamicModelSparse;
      local_temporary_terms = temporary_terms;
    }
  else
    local_output_type = oMatlabDynamicModelSparseLocalTemporaryTerms;

  //----------------------------------------------------------------------
  //For each block
  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {

      //recursive_variables.clear();
      feedback_variables.clear();
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      nze = derivative_endo[block].size();
      nze_other_endo = derivative_other_endo[block].size();
      nze_exo = derivative_exo[block].size();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      unsigned int block_exo_size = exo_block[block].size();
      unsigned int block_exo_det_size = exo_det_block[block].size();
      unsigned int block_other_endo_size = other_endo_block[block].size();
      int block_max_lag = max_leadlag_block[block].first;
      if (global_temporary_terms)
        {
          local_output_type = oMatlabDynamicModelSparse;
          local_temporary_terms = temporary_terms;
        }
      else
        local_output_type = oMatlabDynamicModelSparseLocalTemporaryTerms;

      tmp1_output.str("");
      tmp1_output << dynamic_basename << "_" << block+1 << ".m";
      output.open(tmp1_output.str().c_str(), ios::out | ios::binary);
      output << "%\n";
      output << "% " << tmp1_output.str() << " : Computes dynamic model for Dynare\n";
      output << "%\n";
      output << "% Warning : this file is generated automatically by Dynare\n";
      output << "%           from model file (.mod)\n\n";
      output << "%/\n";
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        {
          output << "function [y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << block+1 << "(y, x, params, jacobian_eval, y_kmin, periods)\n";
        }
      else if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
        output << "function [residual, y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << block+1 << "(y, x, params, it_, jacobian_eval)\n";
      else if (simulation_type == SOLVE_BACKWARD_SIMPLE || simulation_type == SOLVE_FORWARD_SIMPLE)
        output << "function [residual, y, g1, g2, g3, varargout] = " << dynamic_basename << "_" << block+1 << "(y, x, params, it_, jacobian_eval)\n";
      else
        output << "function [residual, y, g1, g2, g3, b, varargout] = " << dynamic_basename << "_" << block+1 << "(y, x, params, periods, jacobian_eval, y_kmin, y_size)\n";
      BlockType block_type;
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        block_type = SIMULTAN;
      else if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE)
        block_type = SIMULTANS;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) < prologue)
        block_type = PROLOGUE;
      else if ((simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_SIMPLE
                || simulation_type == EVALUATE_BACKWARD    || simulation_type == EVALUATE_FORWARD)
               && getBlockFirstEquation(block) >= equations.size() - epilogue)
        block_type = EPILOGUE;
      else
        block_type = SIMULTANS;
      output << "  % ////////////////////////////////////////////////////////////////////////" << endl
             << "  % //" << string("                     Block ").substr(int (log10(block + 1))) << block + 1 << " " << BlockType0(block_type)
             << "          //" << endl
             << "  % //                     Simulation type "
             << BlockSim(simulation_type) << "  //" << endl
             << "  % ////////////////////////////////////////////////////////////////////////" << endl;
      output << "  global options_ oo_;" << endl;
      //The Temporary terms
      if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
        {
          output << "  if(jacobian_eval)\n";
          output << "    g1 = spalloc(" << block_mfs  << ", " << block_mfs*(1+getBlockMaxLag(block)+getBlockMaxLead(block)) << ", " << nze << ");\n";
          output << "    g1_x=spalloc(" << block_size << ", " << (block_exo_size + block_exo_det_size)
            *(1+max(exo_det_max_leadlag_block[block].first, exo_max_leadlag_block[block].first)+max(exo_det_max_leadlag_block[block].second, exo_max_leadlag_block[block].second))
                 << ", " << nze_exo << ");\n";
          output << "    g1_o=spalloc(" << block_size << ", " << block_other_endo_size
            *(1+other_endo_max_leadlag_block[block].first+other_endo_max_leadlag_block[block].second)
                 << ", " << nze_other_endo << ");\n";
          output << "  end;\n";
        }
      else
        {
          output << "  if(jacobian_eval)\n";
          output << "    g1 = spalloc(" << block_size << ", " << block_size*(1+getBlockMaxLag(block)+getBlockMaxLead(block)) << ", " << nze << ");\n";
          output << "    g1_x=spalloc(" << block_size << ", " << (block_exo_size + block_exo_det_size)
            *(1+max(exo_det_max_leadlag_block[block].first, exo_max_leadlag_block[block].first)+max(exo_det_max_leadlag_block[block].second, exo_max_leadlag_block[block].second))
                 << ", " << nze_exo << ");\n";
          output << "    g1_o=spalloc(" << block_size << ", " << block_other_endo_size
            *(1+other_endo_max_leadlag_block[block].first+other_endo_max_leadlag_block[block].second)
                 << ", " << nze_other_endo << ");\n";
          output << "  else\n";
          if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
            {
              output << "    g1 = spalloc(" << block_mfs << "*options_.periods, "
                     << block_mfs << "*(options_.periods+" << max_leadlag_block[block].first+max_leadlag_block[block].second+1 << ")"
                     << ", " << nze << "*options_.periods);\n";
            }
          else
            {
              output << "    g1 = spalloc(" << block_mfs
                     << ", " << block_mfs << ", " << nze << ");\n";
              output << "    g1_tmp_r = spalloc(" << block_recursive
                     << ", " << block_size << ", " << nze << ");\n";
              output << "    g1_tmp_b = spalloc(" << block_mfs
                     << ", " << block_size << ", " << nze << ");\n";
            }
          output << "  end;\n";
        }

      output << "  g2=0;g3=0;\n";
      if (v_temporary_terms_inuse[block].size())
        {
          tmp_output.str("");
          for (temporary_terms_inuse_type::const_iterator it = v_temporary_terms_inuse[block].begin();
               it != v_temporary_terms_inuse[block].end(); it++)
            tmp_output << " T" << *it;
          output << "  global" << tmp_output.str() << ";\n";
        }
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          temporary_terms_type tt2;
          tt2.clear();
          for (int i = 0; i < (int) block_size; i++)
            {
              if (v_temporary_terms[block][i].size() && global_temporary_terms)
                {
                  output << "  " << "% //Temporary variables initialization" << endl
                         << "  " << "T_zeros = zeros(y_kmin+periods, 1);" << endl;
                  for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                       it != v_temporary_terms[block][i].end(); it++)
                    {
                      output << "  ";
                      (*it)->writeOutput(output, oMatlabDynamicModel, local_temporary_terms);
                      output << " = T_zeros;" << endl;
                    }
                }
            }
        }
      if (simulation_type == SOLVE_BACKWARD_SIMPLE || simulation_type == SOLVE_FORWARD_SIMPLE || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        output << "  residual=zeros(" << block_mfs << ",1);\n";
      else if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        output << "  residual=zeros(" << block_mfs << ",y_kmin+periods);\n";
      if (simulation_type == EVALUATE_BACKWARD)
        output << "  for it_ = (y_kmin+periods):y_kmin+1\n";
      if (simulation_type == EVALUATE_FORWARD)
        output << "  for it_ = y_kmin+1:(y_kmin+periods)\n";

      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          output << "  b = zeros(periods*y_size,1);" << endl
                 << "  for it_ = y_kmin+1:(periods+y_kmin)" << endl
                 << "    Per_y_=it_*y_size;" << endl
                 << "    Per_J_=(it_-y_kmin-1)*y_size;" << endl
                 << "    Per_K_=(it_-1)*y_size;" << endl;
          sps = "  ";
        }
      else
        if (simulation_type == EVALUATE_BACKWARD || simulation_type == EVALUATE_FORWARD)
          sps = "  ";
        else
          sps = "";
      // The equations
      for (unsigned int i = 0; i < block_size; i++)
        {
          temporary_terms_type tt2;
          tt2.clear();
          if (v_temporary_terms[block].size())
            {
              output << "  " << "% //Temporary variables" << endl;
              for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                   it != v_temporary_terms[block][i].end(); it++)
                {
                  output << "  " <<  sps;
                  (*it)->writeOutput(output, local_output_type, local_temporary_terms);
                  output << " = ";
                  (*it)->writeOutput(output, local_output_type, tt2);
                  // Insert current node into tt2
                  tt2.insert(*it);
                  output << ";" << endl;
                }
            }

          int variable_ID = getBlockVariableID(block, i);
          int equation_ID = getBlockEquationID(block, i);
          EquationType equ_type = getBlockEquationType(block, i);
          string sModel = symbol_table.getName(symbol_table.getID(eEndogenous, variable_ID));
          eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
          lhs = eq_node->get_arg1();
          rhs = eq_node->get_arg2();
          tmp_output.str("");
          lhs->writeOutput(tmp_output, local_output_type, local_temporary_terms);
          switch (simulation_type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
            evaluation:     if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
                output << "    % equation " << getBlockEquationID(block, i)+1 << " variable : " << sModel
                       << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "    ";
              if (equ_type == E_EVALUATE)
                {
                  output << tmp_output.str();
                  output << " = ";
                  rhs->writeOutput(output, local_output_type, local_temporary_terms);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  output << "%" << tmp_output.str();
                  output << " = ";
                  if (isBlockEquationRenormalized(block, i))
                    {
                      rhs->writeOutput(output, local_output_type, local_temporary_terms);
                      output << "\n    ";
                      tmp_output.str("");
                      eq_node = (BinaryOpNode *) getBlockEquationRenormalizedNodeID(block, i);
                      lhs = eq_node->get_arg1();
                      rhs = eq_node->get_arg2();
                      lhs->writeOutput(output, local_output_type, local_temporary_terms);
                      output << " = ";
                      rhs->writeOutput(output, local_output_type, local_temporary_terms);
                    }
                }
              else
                {
                  cerr << "Type missmatch for equation " << equation_ID+1  << "\n";
                  exit(EXIT_FAILURE);
                }
              output << ";\n";
              break;
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "  % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              output << "  " << "residual(" << i+1-block_recursive << ") = (";
              goto end;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              if (i < block_recursive)
                goto evaluation;
              feedback_variables.push_back(variable_ID);
              output << "    % equation " << equation_ID+1 << " variable : " << sModel
                     << " (" << variable_ID+1 << ") " << c_Equation_Type(equ_type) << endl;
              Uf[equation_ID] << "    b(" << i+1-block_recursive << "+Per_J_) = -residual(" << i+1-block_recursive << ", it_)";
              output << "    residual(" << i+1-block_recursive << ", it_) = (";
              goto end;
            default:
            end:
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, local_output_type, local_temporary_terms);
              output << ");\n";
#ifdef CONDITION
              if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
                output << "  condition(" << i+1 << ")=0;\n";
#endif
            }
        }
      // The Jacobian if we have to solve the block
      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE)
        output << "  " << sps << "% Jacobian  " << endl;
      else
        if (simulation_type == SOLVE_BACKWARD_SIMPLE   || simulation_type == SOLVE_FORWARD_SIMPLE
            || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
          output << "  % Jacobian  " << endl << "  if jacobian_eval" << endl;
        else
          output << "    % Jacobian  " << endl << "    if jacobian_eval" << endl;
      switch (simulation_type)
        {
        case EVALUATE_BACKWARD:
        case EVALUATE_FORWARD:
          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            {
              int lag = it->first.first;
              int eq = it->first.second.first;
              int var = it->first.second.second;
              int eqr = getBlockInitialEquationID(block, eq);
              int varr = getBlockInitialVariableID(block, var);

              NodeID id = it->second;

              output << "      g1(" << eqr+1 << ", " << varr+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            {
              int lag = it->first.first;
              int eq = it->first.second.first;
              int var = it->first.second.second;
              int eqr = getBlockInitialEquationID(block, eq);
              NodeID id = it->second;

              output << "      g1_o(" << eqr+1 << ", " << var+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          output << "      varargout{1}=g1_x;\n";
          output << "      varargout{2}=g1_o;\n";
          output << "    end;" << endl;
          output << "  end;" << endl;
          break;
        case SOLVE_BACKWARD_SIMPLE:
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_COMPLETE:
        case SOLVE_FORWARD_COMPLETE:
          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            {
              int lag = it->first.first;
              unsigned int eq = it->first.second.first;
              unsigned int var = it->first.second.second;
              NodeID id = it->second;

              output << "    g1(" << eq+1 << ", " << var+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }

          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            {
              int lag = it->first.first;
              unsigned int eq = it->first.second.first;
              unsigned int var = it->first.second.second;
              NodeID id = it->second;

              output << "    g1_o(" << eq+1 << ", " << var+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          output << "    varargout{1}=g1_x;\n";
          output << "    varargout{2}=g1_o;\n";
          output << "  else" << endl;
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              unsigned int eq = it->first.first;
              unsigned int var = it->first.second;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              NodeID id = it->second.second;
              int lag = it->second.first;
              if (lag == 0)
                {
                  output << "    g1(" << eq+1 << ", " << var+1-block_recursive << ") = ";
                  id->writeOutput(output, local_output_type, local_temporary_terms);
                  output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                         << "(" << lag
                         << ") " << varr+1
                         << ", equation=" << eqr+1 << endl;
                }

            }
          output << "  end;\n";
          break;
        case SOLVE_TWO_BOUNDARIES_SIMPLE:
        case SOLVE_TWO_BOUNDARIES_COMPLETE:
          output << "    if ~jacobian_eval" << endl;
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
            {
              unsigned int eq = it->first.first;
              unsigned int var = it->first.second;
              unsigned int eqr = getBlockEquationID(block, eq);
              unsigned int varr = getBlockVariableID(block, var);
              ostringstream tmp_output;
              NodeID id = it->second.second;
              int lag = it->second.first;
              if (eq >= block_recursive and var >= block_recursive)
                {
                  if (lag == 0)
                    Uf[eqr] << "+g1(" << eq+1-block_recursive
                            << "+Per_J_, " << var+1-block_recursive
                            << "+Per_K_)*y(it_, " << varr+1 << ")";
                  else if (lag == 1)
                    Uf[eqr] << "+g1(" << eq+1-block_recursive
                            << "+Per_J_, " << var+1-block_recursive
                            << "+Per_y_)*y(it_+1, " << varr+1 << ")";
                  else if (lag > 0)
                    Uf[eqr] << "+g1(" << eq+1-block_recursive
                            << "+Per_J_, " << var+1-block_recursive
                            << "+y_size*(it_+" << lag-1 << "))*y(it_+" << lag << ", " << varr+1 << ")";
                  else if (lag < 0)
                    Uf[eqr] << "+g1(" << eq+1-block_recursive
                            << "+Per_J_, " << var+1-block_recursive
                            << "+y_size*(it_" << lag-1 << "))*y(it_" << lag << ", " << varr+1 << ")";
                  if (lag == 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+Per_K_) = ";
                  else if (lag == 1)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+Per_y_) = ";
                  else if (lag > 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+y_size*(it_+" << lag-1 << ")) = ";
                  else if (lag < 0)
                    tmp_output << "     g1(" << eq+1-block_recursive << "+Per_J_, "
                               << var+1-block_recursive << "+y_size*(it_" << lag-1 << ")) = ";
                  output << " " << tmp_output.str();
                  id->writeOutput(output, local_output_type, local_temporary_terms);
                  output << ";";
                  output << " %2 variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, varr))
                         << "(" << lag << ") " << varr+1
                         << ", equation=" << eqr+1 << " (" << eq+1 << ")" << endl;
                }

#ifdef CONDITION
              output << "  if (fabs(condition[" << eqr << "])<fabs(u[" << u << "+Per_u_]))\n";
              output << "    condition(" << eqr << ")=u(" << u << "+Per_u_);\n";
#endif
            }
          for (unsigned int i = 0; i < block_size; i++)
            {
              if (i >= block_recursive)
                output << "  " << Uf[getBlockEquationID(block, i)].str() << ";\n";
#ifdef CONDITION
              output << "  if (fabs(condition(" << i+1 << "))<fabs(u(" << i << "+Per_u_)))\n";
              output << "    condition(" << i+1 << ")=u(" << i+1 << "+Per_u_);\n";
#endif
            }
#ifdef CONDITION
          for (m = 0; m <= ModelBlock->Block_List[block].Max_Lead+ModelBlock->Block_List[block].Max_Lag; m++)
            {
              k = m-ModelBlock->Block_List[block].Max_Lag;
              for (i = 0; i < ModelBlock->Block_List[block].IM_lead_lag[m].size; i++)
                {
                  unsigned int eq = ModelBlock->Block_List[block].IM_lead_lag[m].Equ_Index[i];
                  unsigned int var = ModelBlock->Block_List[block].IM_lead_lag[m].Var_Index[i];
                  unsigned int u = ModelBlock->Block_List[block].IM_lead_lag[m].u[i];
                  unsigned int eqr = ModelBlock->Block_List[block].IM_lead_lag[m].Equ[i];
                  output << "  u(" << u+1 << "+Per_u_) = u(" << u+1 << "+Per_u_) / condition(" << eqr+1 << ");\n";
                }
            }
          for (i = 0; i < ModelBlock->Block_List[block].Size; i++)
            output << "  u(" << i+1 << "+Per_u_) = u(" << i+1 << "+Per_u_) / condition(" << i+1 << ");\n";
#endif

          output << "    else" << endl;

          for (t_derivative::const_iterator it = derivative_endo[block].begin(); it != derivative_endo[block].end(); it++)
            {
              int lag = it->first.first;
              unsigned int eq = it->first.second.first;
              unsigned int var = it->first.second.second;
              NodeID id = it->second;
              output << "      g1(" << eq+1 << ", " << var+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          for (t_derivative::const_iterator it = derivative_other_endo[block].begin(); it != derivative_other_endo[block].end(); it++)
            {
              int lag = it->first.first;
              unsigned int eq = it->first.second.first;
              unsigned int var = it->first.second.second;
              NodeID id = it->second;

              output << "      g1_o(" << eq+1 << ", " << var+1+(lag+block_max_lag)*block_size << ") = ";
              id->writeOutput(output, local_output_type, local_temporary_terms);
              output << "; % variable=" << symbol_table.getName(symbol_table.getID(eEndogenous, var))
                     << "(" << lag
                     << ") " << var+1
                     << ", equation=" << eq+1 << endl;
            }
          output << "      varargout{1}=g1_x;\n";
          output << "      varargout{2}=g1_o;\n";
          output << "    end;\n";
          output << "  end;\n";
          break;
        default:
          break;
        }
      output.close();
    }
}

void
DynamicModel::writeModelEquationsCodeOrdered(const string file_name, const string bin_basename, map_idx_type map_idx) const
{
  struct Uff_l
  {
    int u, var, lag;
    Uff_l *pNext;
  };

  struct Uff
  {
    Uff_l *Ufl, *Ufl_First;
  };

  int i, v;
  string tmp_s;
  ostringstream tmp_output;
  ofstream code_file;
  NodeID lhs = NULL, rhs = NULL;
  BinaryOpNode *eq_node;
  Uff Uf[symbol_table.endo_nbr()];
  map<NodeID, int> reference_count;
  vector<int> feedback_variables;
  bool file_open = false;

  string main_name = file_name;
  main_name += ".cod";
  code_file.open(main_name.c_str(), ios::out | ios::binary | ios::ate);
  if (!code_file.is_open())
    {
      cout << "Error : Can't open file \"" << main_name << "\" for writing\n";
      exit(EXIT_FAILURE);
    }
  //Temporary variables declaration

  FDIMT_ fdimt(temporary_terms.size());
  fdimt.write(code_file);

  for (unsigned int block = 0; block < getNbBlocks(); block++)
    {
      feedback_variables.clear();
      if (block > 0)
        {
          FENDBLOCK_ fendblock;
          fendblock.write(code_file);
        }
      int count_u;
      int u_count_int = 0;
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      int block_max_lag = max_leadlag_block[block].first;

      if (simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE
          || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          Write_Inf_To_Bin_File(file_name, bin_basename, block, u_count_int, file_open,
                                simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE);
          file_open = true;
        }
      FBEGINBLOCK_ fbeginblock(block_mfs,
                               simulation_type,
                               getBlockFirstEquation(block),
                               block_size,
                               variable_reordered,
                               equation_reordered,
                               blocks_linear[block],
                               symbol_table.endo_nbr(),
                               block_max_lag,
                               block_max_lag,
                               u_count_int
                               );
      fbeginblock.write(code_file);

      // The equations
      for (i = 0; i < (int) block_size; i++)
        {
          //The Temporary terms
          temporary_terms_type tt2;
          tt2.clear();
          if (v_temporary_terms[block][i].size())
            {
              for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
                   it != v_temporary_terms[block][i].end(); it++)
                {
                  (*it)->compile(code_file, false, tt2, map_idx, true, false);
                  FSTPT_ fstpt((int)(map_idx.find((*it)->idx)->second));
                  fstpt.write(code_file);
                  // Insert current node into tt2
                  tt2.insert(*it);
#ifdef DEBUGC
                  cout << "FSTPT " << v << "\n";
                  code_file.write(&FOK, sizeof(FOK));
                  code_file.write(reinterpret_cast<char *>(&k), sizeof(k));
                  ki++;
#endif

                }
            }
#ifdef DEBUGC
          for (temporary_terms_type::const_iterator it = v_temporary_terms[block][i].begin();
               it != v_temporary_terms[block][i].end(); it++)
            {
              map_idx_type::const_iterator ii = map_idx.find((*it)->idx);
              cout << "map_idx[" << (*it)->idx <<"]=" << ii->second << "\n";
            }
#endif

          int variable_ID, equation_ID;
          EquationType equ_type;

          switch (simulation_type)
            {
            evaluation:
            case EVALUATE_BACKWARD:
            case EVALUATE_FORWARD:
              equ_type = getBlockEquationType(block, i);
              if (equ_type == E_EVALUATE)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, true, temporary_terms, map_idx, true, false);
                }
              else if (equ_type == E_EVALUATE_S)
                {
                  eq_node = (BinaryOpNode *) getBlockEquationRenormalizedNodeID(block, i);
                  lhs = eq_node->get_arg1();
                  rhs = eq_node->get_arg2();
                  rhs->compile(code_file, false, temporary_terms, map_idx, true, false);
                  lhs->compile(code_file, true, temporary_terms, map_idx, true, false);
                }
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              if (i < (int) block_recursive)
                goto evaluation;
              variable_ID = getBlockVariableID(block, i);
              equation_ID = getBlockEquationID(block, i);
              feedback_variables.push_back(variable_ID);
              Uf[equation_ID].Ufl = NULL;
              goto end;
            default:
            end:
              eq_node = (BinaryOpNode *) getBlockEquationNodeID(block, i);
              lhs = eq_node->get_arg1();
              rhs = eq_node->get_arg2();
              lhs->compile(code_file, false, temporary_terms, map_idx, true, false);
              rhs->compile(code_file, false, temporary_terms, map_idx, true, false);

              FBINARY_ fbinary(oMinus);
              fbinary.write(code_file);
              FSTPR_ fstpr(i - block_recursive);
              fstpr.write(code_file);
            }
        }
      FENDEQU_ fendequ;
      fendequ.write(code_file);
      // The Jacobian if we have to solve the block
      if    (simulation_type != EVALUATE_BACKWARD
             && simulation_type != EVALUATE_FORWARD)
        {
          switch (simulation_type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FORWARD_SIMPLE:
              compileDerivative(code_file, getBlockEquationID(block, 0), getBlockVariableID(block, 0), 0, map_idx);
              {
                FSTPG_ fstpg(0);
                fstpg.write(code_file);
              }
              break;

            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FORWARD_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
            case SOLVE_TWO_BOUNDARIES_SIMPLE:
              count_u = feedback_variables.size();
              for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != (blocks_derivatives[block]).end(); it++)
                {
                  unsigned int eq = it->first.first;
                  unsigned int var = it->first.second;
                  unsigned int eqr = getBlockEquationID(block, eq);
                  unsigned int varr = getBlockVariableID(block, var);
                  int lag = it->second.first;
                  if (eq >= block_recursive and var >= block_recursive)
                    {
                      if (lag != 0 && (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE))
                        continue;
                      if (!Uf[eqr].Ufl)
                        {
                          Uf[eqr].Ufl = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl_First = Uf[eqr].Ufl;
                        }
                      else
                        {
                          Uf[eqr].Ufl->pNext = (Uff_l *) malloc(sizeof(Uff_l));
                          Uf[eqr].Ufl = Uf[eqr].Ufl->pNext;
                        }
                      Uf[eqr].Ufl->pNext = NULL;
                      Uf[eqr].Ufl->u = count_u;
                      Uf[eqr].Ufl->var = varr;
                      Uf[eqr].Ufl->lag = lag;
                      compileChainRuleDerivative(code_file, eqr, varr, lag, map_idx);
                      FSTPU_ fstpu(count_u);
                      fstpu.write(code_file);
                      count_u++;
                    }
                }
              for (i = 0; i < (int) block_size; i++)
                {
                  if (i >= (int) block_recursive)
                    {
                      FLDR_ fldr(i-block_recursive);
                      fldr.write(code_file);

                      FLDZ_ fldz;
                      fldz.write(code_file);

                      v = getBlockEquationID(block, i);
                      for (Uf[v].Ufl = Uf[v].Ufl_First; Uf[v].Ufl; Uf[v].Ufl = Uf[v].Ufl->pNext)
                        {
                          FLDU_ fldu(Uf[v].Ufl->u);
                          fldu.write(code_file);
                          FLDV_ fldv(eEndogenous, Uf[v].Ufl->var, Uf[v].Ufl->lag);
                          fldv.write(code_file);

                          FBINARY_ fbinary(oTimes);
                          fbinary.write(code_file);

                          FCUML_ fcuml;
                          fcuml.write(code_file);
                        }
                      Uf[v].Ufl = Uf[v].Ufl_First;
                      while (Uf[v].Ufl)
                        {
                          Uf[v].Ufl_First = Uf[v].Ufl->pNext;
                          free(Uf[v].Ufl);
                          Uf[v].Ufl = Uf[v].Ufl_First;
                        }
                      FBINARY_ fbinary(oMinus);
                      fbinary.write(code_file);

                      FSTPU_ fstpu(i - block_recursive);
                      fstpu.write(code_file);
                    }
                }
              break;
            default:
              break;
            }
        }
    }
  FENDBLOCK_ fendblock;
  fendblock.write(code_file);
  FEND_ fend;
  fend.write(code_file);
  code_file.close();
}

void
DynamicModel::writeDynamicMFile(const string &dynamic_basename) const
{
  string filename = dynamic_basename + ".m";

  ofstream mDynamicModelFile;
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "function [residual, g1, g2, g3] = " << dynamic_basename << "(y, x, params, it_)" << endl
                    << "%" << endl
                    << "% Status : Computes dynamic model for Dynare" << endl
                    << "%" << endl
                    << "% Warning : this file is generated automatically by Dynare" << endl
                    << "%           from model file (.mod)" << endl << endl;

  if (isUnaryOpUsed(oSteadyState))
    mDynamicModelFile << "global oo_;" << endl << endl;

  writeDynamicModel(mDynamicModelFile, false);

  mDynamicModelFile.close();
}

void
DynamicModel::writeDynamicCFile(const string &dynamic_basename, const int order) const
{
  string filename = dynamic_basename + ".c";
  ofstream mDynamicModelFile;

  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "/*" << endl
                    << " * " << filename << " : Computes dynamic model for Dynare" << endl
                    << " *" << endl
                    << " * Warning : this file is generated automatically by Dynare" << endl
                    << " *           from model file (.mod)" << endl
                    << endl
                    << " */" << endl
                    << "#include <math.h>" << endl
                    << "#include \"mex.h\"" << endl
                    << endl
                    << "#define max(a, b) (((a) > (b)) ? (a) : (b))" << endl
                    << "#define min(a, b) (((a) > (b)) ? (b) : (a))" << endl;

  // Writing the function body
  writeDynamicModel(mDynamicModelFile, true);

  // Writing the gateway routine
  mDynamicModelFile << "/* The gateway routine */" << endl
                    << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])" << endl
                    << "{" << endl
                    << "  double *y, *x, *params, *steady_state;" << endl
                    << "  double *residual, *g1, *v2, *v3;" << endl
                    << "  int nb_row_x, it_;" << endl
                    << endl
                    << "  /* Check that no derivatives of higher order than computed are being requested */ " << endl
                    << "  if (nlhs > " << order + 1 << ") " << endl
                    << "    mexErrMsgTxt(\"Derivatives of higher order than computed have been requested\"); " << endl
                    << "  /* Create a pointer to the input matrix y. */" << endl
                    << "  y = mxGetPr(prhs[0]);" << endl
                    << endl
                    << "  /* Create a pointer to the input matrix x. */" << endl
                    << "  x = mxGetPr(prhs[1]);" << endl
                    << endl
                    << "  /* Create a pointer to the input matrix params. */" << endl
                    << "  params = mxGetPr(prhs[2]);" << endl
                    << endl
                    << "  /* Fetch time index */" << endl
                    << "  it_ = (int) mxGetScalar(prhs[3]) - 1;" << endl
                    << endl
                    << "  /* Gets number of rows of matrix x. */" << endl
                    << "  nb_row_x = mxGetM(prhs[1]);" << endl
                    << endl
                    << "  residual = NULL;" << endl
                    << "  if (nlhs >= 1)" << endl
                    << "  {" << endl
                    << "     /* Set the output pointer to the output matrix residual. */" << endl
                    << "     plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);" << endl
                    << "     /* Create a C pointer to a copy of the output matrix residual. */" << endl
                    << "     residual = mxGetPr(plhs[0]);" << endl
                    << "  }" << endl
                    << endl
                    << "  g1 = NULL;" << endl
                    << "  if (nlhs >= 2)" << endl
                    << "  {" << endl
                    << "     /* Set the output pointer to the output matrix g1. */" << endl

                    << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << dynJacobianColsNbr << ", mxREAL);" << endl
                    << "     /* Create a C pointer to a copy of the output matrix g1. */" << endl
                    << "     g1 = mxGetPr(plhs[1]);" << endl
                    << "  }" << endl
                    << endl
                    << "  v2 = NULL;" << endl
                    << " if (nlhs >= 3)" << endl
                    << "  {" << endl
                    << "     /* Set the output pointer to the output matrix v2. */" << endl
                    << "     plhs[2] = mxCreateDoubleMatrix(" << NNZDerivatives[1] << ", " << 3
                    << ", mxREAL);" << endl
                    << "     v2 = mxGetPr(plhs[2]);" << endl
                    << "  }" << endl
                    << endl
                    << "  v3 = NULL;" << endl
                    << " if (nlhs >= 4)" << endl
                    << "  {" << endl
                    << "     /* Set the output pointer to the output matrix v3. */" << endl
                    << "     plhs[3] = mxCreateDoubleMatrix(" << NNZDerivatives[2] << ", " << 3 << ", mxREAL);" << endl
                    << "     v3 = mxGetPr(plhs[3]);" << endl
                    << "  }" << endl
                    << endl
                    << "  steady_state = mxGetPr(mxGetField(mexGetVariable(\"global\", \"oo_\"), 0, \"steady_state\"));" << endl
                    << "  /* Call the C subroutines. */" << endl
                    << "  Dynamic(y, x, nb_row_x, params, steady_state, it_, residual, g1, v2, v3);" << endl
                    << "}" << endl;
  mDynamicModelFile.close();
}

string
DynamicModel::reform(const string name1) const
{
  string name = name1;
  int pos = name.find("\\", 0);
  while (pos >= 0)
    {
      if (name.substr(pos + 1, 1) != "\\")
        {
          name = name.insert(pos, "\\");
          pos++;
        }
      pos++;
      pos = name.find("\\", pos);
    }
  return (name);
}

void
DynamicModel::Write_Inf_To_Bin_File(const string &dynamic_basename, const string &bin_basename, const int &num,
                                    int &u_count_int, bool &file_open, bool is_two_boundaries) const
{
  int j;
  std::ofstream SaveCode;
  if (file_open)
    SaveCode.open((bin_basename + "_dynamic.bin").c_str(), ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open((bin_basename + "_dynamic.bin").c_str(), ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << bin_basename << "_dynamic.bin\" for writing\n";
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  unsigned int block_size = getBlockSize(num);
  unsigned int block_mfs = getBlockMfs(num);
  unsigned int block_recursive = block_size - block_mfs;
  for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[num].begin(); it != (blocks_derivatives[num]).end(); it++)
    {
      unsigned int eq = it->first.first;
      unsigned int var = it->first.second;
      int lag = it->second.first;
      if (lag != 0 && !is_two_boundaries)
        continue;
      if (eq >= block_recursive and var >= block_recursive)
        {
          int v = eq - block_recursive;
          SaveCode.write(reinterpret_cast<char *>(&v), sizeof(v));
          int varr = var - block_recursive + lag * block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }

  if (is_two_boundaries)
    u_count_int += block_mfs;
  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int varr = getBlockVariableID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
    }
  for (j = block_recursive; j < (int) block_size; j++)
    {
      unsigned int eqr = getBlockEquationID(num, j);
      SaveCode.write(reinterpret_cast<char *>(&eqr), sizeof(eqr));
    }
  SaveCode.close();
}

void
DynamicModel::writeSparseDynamicMFile(const string &dynamic_basename, const string &basename) const
{
  string sp;
  ofstream mDynamicModelFile;
  ostringstream tmp, tmp1, tmp_eq;
  int prev_Simulation_Type;
  bool OK;
  chdir(basename.c_str());
  string filename = dynamic_basename + ".m";
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  mDynamicModelFile << "%\n";
  mDynamicModelFile << "% " << filename << " : Computes dynamic model for Dynare\n";
  mDynamicModelFile << "%\n";
  mDynamicModelFile << "% Warning : this file is generated automatically by Dynare\n";
  mDynamicModelFile << "%           from model file (.mod)\n\n";
  mDynamicModelFile << "%/\n";

  int Nb_SGE = 0;
  bool skip_head, open_par = false;

  mDynamicModelFile << "function [varargout] = " << dynamic_basename << "(varargin)\n";
  mDynamicModelFile << "  global oo_ options_ M_ ;\n";
  mDynamicModelFile << "  g2=[];g3=[];\n";
  //Temporary variables declaration
  OK = true;
  ostringstream tmp_output;
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    {
      if (OK)
        OK = false;
      else
        tmp_output << " ";
      (*it)->writeOutput(tmp_output, oMatlabStaticModelSparse, temporary_terms);
    }
  if (tmp_output.str().length() > 0)
    mDynamicModelFile << "  global " << tmp_output.str() << " M_ ;\n";

  mDynamicModelFile << "  T_init=zeros(1,options_.periods+M_.maximum_lag+M_.maximum_lead);\n";
  tmp_output.str("");
  for (temporary_terms_type::const_iterator it = temporary_terms.begin();
       it != temporary_terms.end(); it++)
    {
      tmp_output << "  ";
      (*it)->writeOutput(tmp_output, oMatlabDynamicModel, temporary_terms);
      tmp_output << "=T_init;\n";
    }
  if (tmp_output.str().length() > 0)
    mDynamicModelFile << tmp_output.str();

  mDynamicModelFile << "  y_kmin=M_.maximum_lag;" << endl
                    << "  y_kmax=M_.maximum_lead;" << endl
                    << "  y_size=M_.endo_nbr;" << endl
                    << "  if(length(varargin)>0)" << endl
                    << "    %it is a simple evaluation of the dynamic model for time _it" << endl
                    << "    params=varargin{3};" << endl
                    << "    it_=varargin{4};" << endl
                    << "    Per_u_=0;" << endl
                    << "    Per_y_=it_*y_size;" << endl
                    << "    y=varargin{1};" << endl
                    << "    ys=y(it_,:);" << endl
                    << "    x=varargin{2};" << endl;
  prev_Simulation_Type = -1;
  tmp.str("");
  tmp_eq.str("");
  unsigned int nb_blocks = getNbBlocks();
  unsigned int block = 0;
  for (int count_call = 1; block < nb_blocks; block++, count_call++)
    {
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      BlockSimulationType simulation_type = getBlockSimulationType(block);

      if (simulation_type == EVALUATE_FORWARD || simulation_type == EVALUATE_BACKWARD)
        {
          for (unsigned int ik = 0; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
              tmp_eq << " " << getBlockEquationID(block, ik)+1;
            }
        }
      else
        {
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
              tmp_eq << " " << getBlockEquationID(block, ik)+1;
            }
        }
      mDynamicModelFile << "    y_index_eq=[" << tmp_eq.str() << "];\n";
      mDynamicModelFile << "    y_index=[" << tmp.str() << "];\n";

      switch (simulation_type)
        {
        case EVALUATE_FORWARD:
        case EVALUATE_BACKWARD:
          mDynamicModelFile << "    [y, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << block + 1 << "(y, x, params, 1, it_-1, 1);\n";
          mDynamicModelFile << "    residual(y_index_eq)=ys(y_index)-y(it_, y_index);\n";
          break;
        case SOLVE_FORWARD_SIMPLE:
        case SOLVE_BACKWARD_SIMPLE:
          mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << block + 1 << "(y, x, params, it_, 1);\n";
          mDynamicModelFile << "    residual(y_index_eq)=r;\n";
          break;
        case SOLVE_FORWARD_COMPLETE:
        case SOLVE_BACKWARD_COMPLETE:
          mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" << block + 1 << "(y, x, params, it_, 1);\n";
          mDynamicModelFile << "    residual(y_index_eq)=r;\n";
          break;
        case SOLVE_TWO_BOUNDARIES_COMPLETE:
        case SOLVE_TWO_BOUNDARIES_SIMPLE:
          mDynamicModelFile << "    [r, dr(" << count_call << ").g1, dr(" << count_call << ").g2, dr(" << count_call << ").g3, b, dr(" << count_call << ").g1_x, dr(" << count_call << ").g1_o]=" << dynamic_basename << "_" <<  block + 1 << "(y, x, params, it_-" << max_lag << ", 1, " << max_lag << ", " << block_recursive << ");\n";
          mDynamicModelFile << "    residual(y_index_eq)=r(:,M_.maximum_lag+1);\n";
          break;
        default:
          break;
        }
      tmp_eq.str("");
      tmp.str("");
    }
  if (tmp1.str().length())
    {
      mDynamicModelFile << tmp1.str();
      tmp1.str("");
    }
  mDynamicModelFile << "    varargout{1}=residual;" << endl
                    << "    varargout{2}=dr;" << endl
                    << "    return;" << endl
                    << "  end;" << endl
                    << "  %it is the deterministic simulation of the block decomposed dynamic model" << endl
                    << "  if(options_.stack_solve_algo==1)" << endl
                    << "    mthd='Sparse LU';" << endl
                    << "  elseif(options_.stack_solve_algo==2)" << endl
                    << "    mthd='GMRES';" << endl
                    << "  elseif(options_.stack_solve_algo==3)" << endl
                    << "    mthd='BICGSTAB';" << endl
                    << "  elseif(options_.stack_solve_algo==4)" << endl
                    << "    mthd='OPTIMPATH';" << endl
                    << "  else" << endl
                    << "    mthd='UNKNOWN';" << endl
                    << "  end;" << endl
                    << "  disp (['-----------------------------------------------------']) ;" << endl
                    << "  disp (['MODEL SIMULATION: (method=' mthd ')']) ;" << endl
                    << "  fprintf('\\n') ;" << endl
                    << "  periods=options_.periods;" << endl
                    << "  maxit_=options_.maxit_;" << endl
                    << "  solve_tolf=options_.solve_tolf;" << endl
                    << "  y=oo_.endo_simul';" << endl
                    << "  x=oo_.exo_simul;" << endl;

  prev_Simulation_Type = -1;
  mDynamicModelFile << "  params=M_.params;\n";
  mDynamicModelFile << "  oo_.deterministic_simulation.status = 0;\n";
  for (block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      unsigned int block_mfs = getBlockMfs(block);
      unsigned int block_recursive = block_size - block_mfs;
      BlockSimulationType simulation_type = getBlockSimulationType(block);

      if (BlockSim(prev_Simulation_Type) == BlockSim(simulation_type)
          && (simulation_type == EVALUATE_FORWARD || simulation_type == EVALUATE_BACKWARD))
        skip_head = true;
      else
        skip_head = false;
      if ((simulation_type == EVALUATE_FORWARD) && (block_size))
        {
          if (!skip_head)
            {
              if (open_par)
                {
                  mDynamicModelFile << "  end\n";
                }
              mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.error = 0;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.iterations = 0;\n";
              mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
              mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
              mDynamicModelFile << "  else\n";
              mDynamicModelFile << "    blck_num = 1;\n";
              mDynamicModelFile << "  end;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).status = 1;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).error = 0;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;\n";
              mDynamicModelFile << "  g1=[];g2=[];g3=[];\n";
              mDynamicModelFile << "  y=" << dynamic_basename << "_" << block + 1 << "(y, x, params, 0, y_kmin, periods);\n";
              mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);\n";
              mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
              mDynamicModelFile << "    disp(['Inf or Nan value during the evaluation of block " << block <<"']);\n";
              mDynamicModelFile << "    return;\n";
              mDynamicModelFile << "  end;\n";
            }
        }
      else if ((simulation_type == EVALUATE_BACKWARD) && (block_size))
        {
          if (!skip_head)
            {
              if (open_par)
                {
                  mDynamicModelFile << "  end\n";
                }
              mDynamicModelFile << "  oo_.deterministic_simulation.status = 1;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.error = 0;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.iterations = 0;\n";
              mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
              mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
              mDynamicModelFile << "  else\n";
              mDynamicModelFile << "    blck_num = 1;\n";
              mDynamicModelFile << "  end;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).status = 1;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).error = 0;\n";
              mDynamicModelFile << "  oo_.deterministic_simulation.block(blck_num).iterations = 0;\n";
              mDynamicModelFile << "  g1=[];g2=[];g3=[];\n";
              mDynamicModelFile << "  " << dynamic_basename << "_" << block + 1 << "(y, x, params, 0, y_kmin, periods);\n";
              mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);\n";
              mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
              mDynamicModelFile << "    disp(['Inf or Nan value during the evaluation of block " << block <<"']);\n";
              mDynamicModelFile << "    return;\n";
              mDynamicModelFile << "  end;\n";
            }
        }
      else if ((simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_FORWARD_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par = false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          tmp.str("");
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
            }
          mDynamicModelFile << "  y_index = [" << tmp.str() << "];\n";
          int nze = blocks_derivatives[block].size();
          mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
          mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
          mDynamicModelFile << "  else\n";
          mDynamicModelFile << "    blck_num = 1;\n";
          mDynamicModelFile << "  end;\n";
          mDynamicModelFile << "  y = solve_one_boundary('"  << dynamic_basename << "_" <<  block + 1 << "'"
                            <<", y, x, params, y_index, " << nze
                            <<", options_.periods, " << blocks_linear[block]
                            <<", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);\n";
          mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);\n";
          mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
          mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);\n";
          mDynamicModelFile << "    return;\n";
          mDynamicModelFile << "  end;\n";
        }
      else if ((simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par = false;
          mDynamicModelFile << "  g1=0;\n";
          mDynamicModelFile << "  r=0;\n";
          tmp.str("");
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            {
              tmp << " " << getBlockVariableID(block, ik)+1;
            }
          mDynamicModelFile << "  y_index = [" << tmp.str() << "];\n";
          int nze = blocks_derivatives[block].size();

          mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
          mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
          mDynamicModelFile << "  else\n";
          mDynamicModelFile << "    blck_num = 1;\n";
          mDynamicModelFile << "  end;\n";
          mDynamicModelFile << "  y = solve_one_boundary('"  << dynamic_basename << "_" <<  block + 1 << "'"
                            <<", y, x, params, y_index, " << nze
                            <<", options_.periods, " << blocks_linear[block]
                            <<", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo, 1, 1, 0);\n";
          mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);\n";
          mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
          mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);\n";
          mDynamicModelFile << "    return;\n";
          mDynamicModelFile << "  end;\n";
        }
      else if ((simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_size))
        {
          if (open_par)
            mDynamicModelFile << "  end\n";
          open_par = false;
          Nb_SGE++;
          int nze = blocks_derivatives[block].size();
          mDynamicModelFile << "  y_index=[";
          for (unsigned int ik = block_recursive; ik < block_size; ik++)
            {
              mDynamicModelFile << " " << getBlockVariableID(block, ik)+1;
            }
          mDynamicModelFile << "  ];\n";
          mDynamicModelFile << "  if(isfield(oo_.deterministic_simulation,'block'))\n";
          mDynamicModelFile << "    blck_num = length(oo_.deterministic_simulation.block)+1;\n";
          mDynamicModelFile << "  else\n";
          mDynamicModelFile << "    blck_num = 1;\n";
          mDynamicModelFile << "  end;\n";
          mDynamicModelFile << "  y = solve_two_boundaries('" << dynamic_basename << "_" <<  block + 1 << "'"
                            <<", y, x, params, y_index, " << nze
                            <<", options_.periods, " << max_leadlag_block[block].first
                            <<", " << max_leadlag_block[block].second
                            <<", " << blocks_linear[block]
                            <<", blck_num, y_kmin, options_.maxit_, options_.solve_tolf, options_.slowc, " << cutoff << ", options_.stack_solve_algo);\n";
          mDynamicModelFile << "  tmp = y(:,M_.block_structure.block(" << block + 1 << ").variable);\n";
          mDynamicModelFile << "  if(isnan(tmp) | isinf(tmp))\n";
          mDynamicModelFile << "    disp(['Inf or Nan value during the resolution of block " << block <<"']);\n";
          mDynamicModelFile << "    return;\n";
          mDynamicModelFile << "  end;\n";
        }
      prev_Simulation_Type = simulation_type;
    }
  if (open_par)
    mDynamicModelFile << "  end;\n";
  open_par = false;
  mDynamicModelFile << "  oo_.endo_simul = y';\n";
  mDynamicModelFile << "return;\n";

  mDynamicModelFile.close();

  writeModelEquationsOrdered_M(dynamic_basename);

  chdir("..");
}

void
DynamicModel::writeDynamicModel(ostream &DynamicOutput, bool use_dll) const
{
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;  // Used for storing Hessian equations
  ostringstream third_derivatives_output;

  ExprNodeOutputType output_type = (use_dll ? oCDynamicModel : oMatlabDynamicModel);

  writeModelLocalVariables(model_output, output_type);

  writeTemporaryTerms(temporary_terms, model_output, output_type);

  writeModelEquations(model_output, output_type);

  int nrows = equations.size();
  int hessianColsNbr = dynJacobianColsNbr * dynJacobianColsNbr;

  // Writing Jacobian
  for (first_derivatives_type::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var = it->first.second;
      NodeID d1 = it->second;

      jacobian_output << "g1";
      jacobianHelper(jacobian_output, eq, getDynJacobianCol(var), output_type);
      jacobian_output << "=";
      d1->writeOutput(jacobian_output, output_type, temporary_terms);
      jacobian_output << ";" << endl;
    }

  // Writing Hessian
  int k = 0; // Keep the line of a 2nd derivative in v2
  for (second_derivatives_type::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second.first;
      int var2 = it->first.second.second;
      NodeID d2 = it->second;

      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);

      int col_nb = id1 * dynJacobianColsNbr + id2;
      int col_nb_sym = id2 * dynJacobianColsNbr + id1;

      sparseHelper(2, hessian_output, k, 0, output_type);
      hessian_output << "=" << eq + 1 << ";" << endl;

      sparseHelper(2, hessian_output, k, 1, output_type);
      hessian_output << "=" << col_nb + 1 << ";" << endl;

      sparseHelper(2, hessian_output, k, 2, output_type);
      hessian_output << "=";
      d2->writeOutput(hessian_output, output_type, temporary_terms);
      hessian_output << ";" << endl;

      k++;

      // Treating symetric elements
      if (id1 != id2)
        {
          sparseHelper(2, hessian_output, k, 0, output_type);
          hessian_output << "=" << eq + 1 << ";" << endl;

          sparseHelper(2, hessian_output, k, 1, output_type);
          hessian_output << "=" << col_nb_sym + 1 << ";" << endl;

          sparseHelper(2, hessian_output, k, 2, output_type);
          hessian_output << "=";
          sparseHelper(2, hessian_output, k-1, 2, output_type);
          hessian_output << ";" << endl;

          k++;
        }
    }

  // Writing third derivatives
  k = 0; // Keep the line of a 3rd derivative in v3
  for (third_derivatives_type::const_iterator it = third_derivatives.begin();
       it != third_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second.first;
      int var2 = it->first.second.second.first;
      int var3 = it->first.second.second.second;
      NodeID d3 = it->second;

      int id1 = getDynJacobianCol(var1);
      int id2 = getDynJacobianCol(var2);
      int id3 = getDynJacobianCol(var3);

      // Reference column number for the g3 matrix
      int ref_col = id1 * hessianColsNbr + id2 * dynJacobianColsNbr + id3;

      sparseHelper(3, third_derivatives_output, k, 0, output_type);
      third_derivatives_output << "=" << eq + 1 << ";" << endl;

      sparseHelper(3, third_derivatives_output, k, 1, output_type);
      third_derivatives_output << "=" << ref_col + 1 << ";" << endl;

      sparseHelper(3, third_derivatives_output, k, 2, output_type);
      third_derivatives_output << "=";
      d3->writeOutput(third_derivatives_output, output_type, temporary_terms);
      third_derivatives_output << ";" << endl;

      k++;

      // Compute the column numbers for the 5 other permutations of (id1,id2,id3) and store them in a set (to avoid duplicates if two indexes are equal)
      set<int> cols;
      cols.insert(id1 * hessianColsNbr + id3 * dynJacobianColsNbr + id2);
      cols.insert(id2 * hessianColsNbr + id1 * dynJacobianColsNbr + id3);
      cols.insert(id2 * hessianColsNbr + id3 * dynJacobianColsNbr + id1);
      cols.insert(id3 * hessianColsNbr + id1 * dynJacobianColsNbr + id2);
      cols.insert(id3 * hessianColsNbr + id2 * dynJacobianColsNbr + id1);

      int k2 = 0; // Keeps the offset of the permutation relative to k
      for (set<int>::iterator it2 = cols.begin(); it2 != cols.end(); it2++)
        if (*it2 != ref_col)
          {
            sparseHelper(3, third_derivatives_output, k+k2, 0, output_type);
            third_derivatives_output << "=" << eq + 1 << ";" << endl;

            sparseHelper(3, third_derivatives_output, k+k2, 1, output_type);
            third_derivatives_output << "=" << *it2 + 1 << ";" << endl;

            sparseHelper(3, third_derivatives_output, k+k2, 2, output_type);
            third_derivatives_output << "=";
            sparseHelper(3, third_derivatives_output, k, 2, output_type);
            third_derivatives_output << ";" << endl;

            k2++;
          }
      k += k2;
    }

  if (!use_dll)
    {
      DynamicOutput << "%" << endl
                    << "% Model equations" << endl
                    << "%" << endl
                    << endl
                    << "residual = zeros(" << nrows << ", 1);" << endl
                    << model_output.str()
        // Writing initialization instruction for matrix g1
                    << "if nargout >= 2," << endl
                    << "  g1 = zeros(" << nrows << ", " << dynJacobianColsNbr << ");" << endl
                    << endl
                    << "%" << endl
                    << "% Jacobian matrix" << endl
                    << "%" << endl
                    << endl
                    << jacobian_output.str()
                    << "end" << endl;

      // Initialize g2 matrix
      DynamicOutput << "if nargout >= 3," << endl
                    << "%" << endl
                    << "% Hessian matrix" << endl
                    << "%" << endl
                    << endl;
      if (second_derivatives.size())
        DynamicOutput << "  v2 = zeros(" << NNZDerivatives[1] << ",3);" << endl
                      << hessian_output.str()
                      << "  g2 = sparse(v2(:,1),v2(:,2),v2(:,3)," << nrows << "," << hessianColsNbr << ");" << endl;
      else // Either hessian is all zero, or we didn't compute it
        DynamicOutput << "  g2 = sparse([],[],[]," << nrows << "," << hessianColsNbr << ");" << endl;
      DynamicOutput << "end;" << endl;

      // Initialize g3 matrix
      DynamicOutput << "if nargout >= 4," << endl
                    << "%" << endl
                    << "% Third order derivatives" << endl
                    << "%" << endl
                    << endl;
      int ncols = hessianColsNbr * dynJacobianColsNbr;
      if (third_derivatives.size())
        DynamicOutput << "  v3 = zeros(" << NNZDerivatives[2] << ",3);" << endl
                      << third_derivatives_output.str()
                      << "  g3 = sparse(v3(:,1),v3(:,2),v3(:,3)," << nrows << "," << ncols << ");" << endl;
      else // Either 3rd derivatives is all zero, or we didn't compute it
        DynamicOutput << "  g3 = sparse([],[],[]," << nrows << "," << ncols << ");" << endl;

      DynamicOutput << "end;" << endl;
    }
  else
    {
      DynamicOutput << "void Dynamic(double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, double *g1, double *v2, double *v3)" << endl
                    << "{" << endl
                    << "  double lhs, rhs;" << endl
                    << endl
                    << "  /* Residual equations */" << endl
                    << model_output.str()
                    << "  /* Jacobian  */" << endl
                    << "  if (g1 == NULL)" << endl
                    << "    return;" << endl
                    << "  else" << endl
                    << "    {" << endl
                    << jacobian_output.str()
                    << "    }" << endl;

      if (second_derivatives.size())
        DynamicOutput << "  /* Hessian for endogenous and exogenous variables */" << endl
                      << "  if (v2 == NULL)" << endl
                      << "    return;" << endl
                      << "  else" << endl
                      << "    {" << endl
                      << hessian_output.str()
                      << "    }" << endl;

      if (third_derivatives.size())
        DynamicOutput << "  /* Third derivatives for endogenous and exogenous variables */" << endl
                      << "  if (v3 == NULL)" << endl
                      << "    return;" << endl
                      << "  else" << endl
                      << "    {" << endl
                      << third_derivatives_output.str()
                      << "    }" << endl;

      DynamicOutput << "}" << endl << endl;
    }
}

void
DynamicModel::writeOutput(ostream &output, const string &basename, bool block_decomposition, bool byte_code, bool use_dll, int order) const
{
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */

  output << "M_.lead_lag_incidence = [";
  // Loop on endogenous variables
  for (int endoID = 0; endoID < symbol_table.endo_nbr(); endoID++)
    {
      output << endl;
      // Loop on periods
      for (int lag = -max_endo_lag; lag <= max_endo_lead; lag++)
        {
          // Print variableID if exists with current period, otherwise print 0
          try
            {
              int varID = getDerivID(symbol_table.getID(eEndogenous, endoID), lag);
              output << " " << getDynJacobianCol(varID) + 1;
            }
          catch (UnknownDerivIDException &e)
            {
              output << " 0";
            }
        }
      output << ";";
    }
  output << "]';" << endl;

  // Write equation tags
  output << "M_.equations_tags = {" << endl;
  for (unsigned int i = 0; i < equation_tags.size(); i++)
    output << "  " << equation_tags[i].first + 1 << " , '"
           << equation_tags[i].second.first << "' , '"
           << equation_tags[i].second.second << "' ;" << endl;
  output << "};" << endl;

  //In case of sparse model, writes the block_decomposition structure of the model
  if (block_decomposition)
    {
      int count_lead_lag_incidence = 0;
      int max_lead, max_lag, max_lag_endo, max_lead_endo, max_lag_exo, max_lead_exo;
      unsigned int nb_blocks = getNbBlocks();
      for (unsigned int block = 0; block < nb_blocks; block++)
        {
          //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
          count_lead_lag_incidence = 0;
          BlockSimulationType simulation_type = getBlockSimulationType(block);
          int block_size = getBlockSize(block);
          max_lag  = max_leadlag_block[block].first;
          max_lead = max_leadlag_block[block].second;
          max_lag_endo = endo_max_leadlag_block[block].first;
          max_lead_endo = endo_max_leadlag_block[block].second;
          max_lag_exo = max(exo_max_leadlag_block[block].first, exo_det_max_leadlag_block[block].first);
          max_lead_exo = max(exo_max_leadlag_block[block].second, exo_det_max_leadlag_block[block].second);
          vector<int> exogenous(symbol_table.exo_nbr(), -1);
          vector<int>::iterator it_exogenous;
          exogenous.clear();
          ostringstream tmp_s, tmp_s_eq;
          tmp_s.str("");
          tmp_s_eq.str("");
          for (int i = 0; i < block_size; i++)
            {
              tmp_s << " " << getBlockVariableID(block, i)+1;
              tmp_s_eq << " " << getBlockEquationID(block, i)+1;
            }
          it_exogenous = exogenous.begin();
          for (t_lag_var::const_iterator it = exo_block[block].begin(); it != exo_block[block].end(); it++)
            it_exogenous = set_union(it->second.begin(), it->second.end(), exogenous.begin(), exogenous.end(), it_exogenous);
          output << "M_.block_structure.block(" << block+1 << ").num = " << block+1 << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").Simulation_Type = " << simulation_type << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_lag = " << max_lag << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_lead = " << max_lead << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_endo_lag = " << max_lag_endo << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_endo_lead = " << max_lead_endo << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_exo_lag = " << max_lag_exo << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").maximum_exo_lead = " << max_lead_exo << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").endo_nbr = " << block_size << ";\n";
          output << "M_.block_structure.block(" << block+1 << ").equation = [" << tmp_s_eq.str() << "];\n";
          output << "M_.block_structure.block(" << block+1 << ").variable = [" << tmp_s.str() << "];\n";
          output << "M_.block_structure.block(" << block+1 << ").exogenous = [";
          int i = 0;
          for (it_exogenous = exogenous.begin(); it_exogenous != exogenous.end(); it_exogenous++)
            if (*it_exogenous >= 0)
              {
                output << " " << *it_exogenous+1;
                i++;
              }
          output << "];\n";
          output << "M_.block_structure.block(" << block+1 << ").exo_nbr = " << i << ";\n";

          output << "M_.block_structure.block(" << block+1 << ").exo_det_nbr = " << exo_det_block.size() << ";\n";

          tmp_s.str("");
          count_lead_lag_incidence = 0;
          dynamic_jacob_map reordered_dynamic_jacobian;
          for (t_block_derivatives_equation_variable_laglead_nodeid::const_iterator it = blocks_derivatives[block].begin(); it != blocks_derivatives[block].end(); it++)
            reordered_dynamic_jacobian[make_pair(it->second.first, make_pair(it->first.second, it->first.first))] = it->second.second;
          output << "M_.block_structure.block(" << block+1 << ").lead_lag_incidence = [];\n";
          int last_var = -1;
          for (int lag = -max_lag_endo; lag < max_lead_endo+1; lag++)
            {
              last_var = 0;
              for (dynamic_jacob_map::const_iterator it = reordered_dynamic_jacobian.begin(); it != reordered_dynamic_jacobian.end(); it++)
                {
                  if (lag == it->first.first && last_var != it->first.second.first)
                    {
                      count_lead_lag_incidence++;
                      for (int i = last_var; i < it->first.second.first-1; i++)
                        tmp_s << " 0";
                      if (tmp_s.str().length())
                        tmp_s << " ";
                      tmp_s << count_lead_lag_incidence;
                      last_var = it->first.second.first;
                    }
                }
              for (int i = last_var + 1; i < block_size; i++)
                tmp_s << " 0";
              output << "M_.block_structure.block(" << block+1 << ").lead_lag_incidence = [ M_.block_structure.block(" << block+1 << ").lead_lag_incidence; " << tmp_s.str() << "]; %lag = " << lag << "\n";
              tmp_s.str("");
            }
        }
      string cst_s;
      int nb_endo = symbol_table.endo_nbr();
      output << "M_.block_structure.variable_reordered = [";
      for (int i = 0; i < nb_endo; i++)
        output << " " << variable_reordered[i]+1;
      output << "];\n";
      output << "M_.block_structure.equation_reordered = [";
      for (int i = 0; i < nb_endo; i++)
        output << " " << equation_reordered[i]+1;
      output << "];\n";
      map<pair< int, pair<int, int> >,  int>  lag_row_incidence;
      for (first_derivatives_type::const_iterator it = first_derivatives.begin();
           it != first_derivatives.end(); it++)
        {
          int deriv_id = it->first.second;
          if (getTypeByDerivID(deriv_id) == eEndogenous)
            {
              int eq = it->first.first;
              int symb = getSymbIDByDerivID(deriv_id);
              int var = symbol_table.getTypeSpecificID(symb);
              int lag = getLagByDerivID(deriv_id);
              int eqr = inv_equation_reordered[eq];
              int varr = inv_variable_reordered[var];
              lag_row_incidence[make_pair(lag, make_pair(eqr, varr))] = 1;
            }
        }
      int prev_lag = -1000000;
      for (map<pair< int, pair<int, int> >,  int>::const_iterator it = lag_row_incidence.begin(); it != lag_row_incidence.end(); it++)
        {
          if (prev_lag != it->first.first)
            {
              if (prev_lag != -1000000)
                output << "];\n";
              prev_lag = it->first.first;
              output << "M_.block_structure.incidence(" << max_endo_lag+it->first.first+1 << ").lead_lag = " << prev_lag << ";\n";
              output << "M_.block_structure.incidence(" << max_endo_lag+it->first.first+1 << ").sparse_IM = [";
            }
          output << it->first.second.first << " " << it->first.second.second << ";\n";
        }
      output << "];\n";
    }
  // Writing initialization for some other variables
  output << "M_.exo_names_orig_ord = [1:" << symbol_table.exo_nbr() << "];" << endl
         << "M_.maximum_lag = " << max_lag << ";" << endl
         << "M_.maximum_lead = " << max_lead << ";" << endl;
  if (symbol_table.endo_nbr())
    {
      output << "M_.maximum_endo_lag = " << max_endo_lag << ";" << endl
             << "M_.maximum_endo_lead = " << max_endo_lead << ";" << endl
             << "oo_.steady_state = zeros(" << symbol_table.endo_nbr() << ", 1);" << endl;
    }
  if (symbol_table.exo_nbr())
    {
      output << "M_.maximum_exo_lag = " << max_exo_lag << ";" << endl
             << "M_.maximum_exo_lead = " << max_exo_lead << ";" << endl
             << "oo_.exo_steady_state = zeros(" << symbol_table.exo_nbr() << ", 1);" << endl;
    }
  if (symbol_table.exo_det_nbr())
    {
      output << "M_.maximum_exo_det_lag = " << max_exo_det_lag << ";" << endl
             << "M_.maximum_exo_det_lead = " << max_exo_det_lead << ";" << endl
             << "oo_.exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr() << ", 1);" << endl;
    }
  if (symbol_table.param_nbr())
    output << "M_.params = repmat(NaN," << symbol_table.param_nbr() << ", 1);" << endl;

  // Write number of non-zero derivatives
  // Use -1 if the derivatives have not been computed
  output << "M_.NNZDerivatives = zeros(3, 1);" << endl
         << "M_.NNZDerivatives(1) = " << NNZDerivatives[0] << ";" << endl;
  if (order > 1)
    {
      output << "M_.NNZDerivatives(2) = " << NNZDerivatives[1] << ";" << endl;
      if (order > 2)
	output << "M_.NNZDerivatives(3) = " << NNZDerivatives[2] << ";" << endl;
      else
	output << "M_.NNZDerivatives(3) = -1;" << endl;
    }
  else
    output << "M_.NNZDerivatives(2) = -1;" << endl
	   << "M_.NNZDerivatives(3) = -1;" << endl;

}

map<pair<int, pair<int, int > >, NodeID>
DynamicModel::collect_first_order_derivatives_endogenous()
{
  map<pair<int, pair<int, int > >, NodeID> endo_derivatives;
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      if (getTypeByDerivID(it2->first.second) == eEndogenous)
        {
          int eq = it2->first.first;
          int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(it2->first.second));
          int lag = getLagByDerivID(it2->first.second);
          endo_derivatives[make_pair(eq, make_pair(var, lag))] = it2->second;
        }
    }
  return endo_derivatives;
}

void
DynamicModel::computingPass(bool jacobianExo, bool hessian, bool thirdDerivatives, bool paramsDerivatives,
                            const eval_context_type &eval_context, bool no_tmp_terms, bool block, bool use_dll)
{
  assert(jacobianExo || !(hessian || thirdDerivatives || paramsDerivatives));

  // Prepare for derivation
  computeDerivIDs();

  // Computes dynamic jacobian columns, must be done after computeDerivIDs()
  computeDynJacobianCols(jacobianExo);

  // Compute derivatives w.r. to all endogenous, and possibly exogenous and exogenous deterministic
  set<int> vars;
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first.first);
      if (type == eEndogenous || (jacobianExo && (type == eExogenous || type == eExogenousDet)))
        vars.insert(it->second);
    }

  // Launch computations
  cout << "Computing dynamic model derivatives:" << endl
       << " - order 1" << endl;
  computeJacobian(vars);

  if (hessian)
    {
      cout << " - order 2" << endl;
      computeHessian(vars);
    }

  if (paramsDerivatives)
    {
      cout << " - order 2 (derivatives of Jacobian w.r. to parameters)" << endl;
      computeParamsDerivatives();

      if (!no_tmp_terms)
        computeParamsDerivativesTemporaryTerms();
    }

  if (thirdDerivatives)
    {
      cout << " - order 3" << endl;
      computeThirdDerivatives(vars);
    }

  if (block)
    {
      jacob_map contemporaneous_jacobian, static_jacobian;

      // for each block contains pair<Size, Feddback_variable>
      vector<pair<int, int> > blocks;

      evaluateAndReduceJacobian(eval_context, contemporaneous_jacobian, static_jacobian, dynamic_jacobian, cutoff, false);

      computeNonSingularNormalization(contemporaneous_jacobian, cutoff, static_jacobian, dynamic_jacobian);

      computePrologueAndEpilogue(static_jacobian, equation_reordered, variable_reordered, prologue, epilogue);

      map<pair<int, pair<int, int> >, NodeID> first_order_endo_derivatives = collect_first_order_derivatives_endogenous();

      equation_type_and_normalized_equation = equationTypeDetermination(equations, first_order_endo_derivatives, variable_reordered, equation_reordered, mfs);

      cout << "Finding the optimal block decomposition of the model ...\n";

      if (prologue+epilogue < (unsigned int) equation_number())
        computeBlockDecompositionAndFeedbackVariablesForEachBlock(static_jacobian, dynamic_jacobian, prologue, epilogue, equation_reordered, variable_reordered, blocks, equation_type_and_normalized_equation, false, true, mfs, inv_equation_reordered, inv_variable_reordered);

      block_type_firstequation_size_mfs = reduceBlocksAndTypeDetermination(dynamic_jacobian, prologue, epilogue, blocks, equations, equation_type_and_normalized_equation, variable_reordered, equation_reordered);

      printBlockDecomposition(blocks);

      computeChainRuleJacobian(blocks_derivatives);

      blocks_linear = BlockLinear(blocks_derivatives, variable_reordered);

      collect_block_first_order_derivatives();

      global_temporary_terms = true;
      if (!no_tmp_terms)
        computeTemporaryTermsOrdered();

    }
  else
    if (!no_tmp_terms)
      computeTemporaryTerms(!use_dll);
}

map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>
DynamicModel::get_Derivatives(int block)
{
  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives;
  Derivatives.clear();
  int max_lag  = getBlockMaxLag(block);
  int max_lead = getBlockMaxLead(block);
  int block_size = getBlockSize(block);
  int block_nb_recursive = block_size - getBlockMfs(block);
  for (int lag = -max_lag; lag <= max_lead; lag++)
    {
      for (int eq = 0; eq < block_size; eq++)
        {
          int eqr = getBlockEquationID(block, eq);
          for (int var = 0; var < block_size; var++)
            {
              int varr = getBlockVariableID(block, var);
              if (dynamic_jacobian.find(make_pair(lag, make_pair(eqr, varr))) != dynamic_jacobian.end())
                {
                  bool OK = true;
                  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator its = Derivatives.find(make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr)));
                  if (its != Derivatives.end())
                    {
                      if (its->second == 2)
                        OK = false;
                    }

                  if (OK)
                    {
                      if (getBlockEquationType(block, eq) == E_EVALUATE_S and eq < block_nb_recursive)
                        //It's a normalized equation, we have to recompute the derivative using chain rule derivative function
                        Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 1;
                      else
                        //It's a feedback equation we can use the derivatives
                        Derivatives[make_pair(make_pair(lag, make_pair(eq, var)), make_pair(eqr, varr))] = 0;
                    }
                  if (var < block_nb_recursive)
                    {
                      int eqs = getBlockEquationID(block, var);
                      for (int vars = block_nb_recursive; vars < block_size; vars++)
                        {
                          int varrs = getBlockVariableID(block, vars);
                          //A new derivative needs to be computed using the chain rule derivative function (a feedback variable appears in a recursive equation)
                          if (Derivatives.find(make_pair(make_pair(lag, make_pair(var, vars)), make_pair(eqs, varrs))) != Derivatives.end())
                            Derivatives[make_pair(make_pair(lag, make_pair(eq, vars)), make_pair(eqr, varrs))] = 2;
                        }
                    }
                }
            }
        }
    }
  return (Derivatives);
}

void
DynamicModel::computeChainRuleJacobian(t_blocks_derivatives &blocks_derivatives)
{
  map<int, NodeID> recursive_variables;
  unsigned int nb_blocks = getNbBlocks();
  blocks_derivatives = t_blocks_derivatives(nb_blocks);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      t_block_derivatives_equation_variable_laglead_nodeid tmp_derivatives;
      recursive_variables.clear();
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      int block_nb_mfs = getBlockMfs(block);
      int block_nb_recursives = block_size - block_nb_mfs;
      if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE or simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          blocks_derivatives.push_back(t_block_derivatives_equation_variable_laglead_nodeid(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedNodeID(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationNodeID(block, i);
            }
          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> Derivatives = get_Derivatives(block);
          map<pair<pair<int, pair<int, int> >, pair<int, int> >, int>::const_iterator it = Derivatives.begin();
          for (int i = 0; i < (int) Derivatives.size(); i++)
            {
              int Deriv_type = it->second;
              pair<pair<int, pair<int, int> >, pair<int, int> > it_l(it->first);
              it++;
              int lag = it_l.first.first;
              int eq = it_l.first.second.first;
              int var = it_l.first.second.second;
              int eqr = it_l.second.first;
              int varr = it_l.second.second;
              if (Deriv_type == 0)
                first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = first_derivatives[make_pair(eqr, getDerivID(symbol_table.getID(eEndogenous, varr), lag))];
              else if (Deriv_type == 1)
                first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
              else if (Deriv_type == 2)
                {
                  if (getBlockEquationType(block, eq) == E_EVALUATE_S && eq < block_nb_recursives)
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = (equation_type_and_normalized_equation[eqr].second)->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                  else
                    first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))] = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), lag), recursive_variables);
                }
              tmp_derivatives.push_back(make_pair(make_pair(eq, var), make_pair(lag, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, lag))])));
            }
        }
      else if (simulation_type == SOLVE_BACKWARD_SIMPLE or simulation_type == SOLVE_FORWARD_SIMPLE
               or simulation_type == SOLVE_BACKWARD_COMPLETE or simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          blocks_derivatives.push_back(t_block_derivatives_equation_variable_laglead_nodeid(0));
          for (int i = 0; i < block_nb_recursives; i++)
            {
              if (getBlockEquationType(block, i) == E_EVALUATE_S)
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationRenormalizedNodeID(block, i);
              else
                recursive_variables[getDerivID(symbol_table.getID(eEndogenous, getBlockVariableID(block, i)), 0)] = getBlockEquationNodeID(block, i);
            }
          for (int eq = block_nb_recursives; eq < block_size; eq++)
            {
              int eqr = getBlockEquationID(block, eq);
              for (int var = block_nb_recursives; var < block_size; var++)
                {
                  int varr = getBlockVariableID(block, var);
                  NodeID d1 = equations[eqr]->getChainRuleDerivative(getDerivID(symbol_table.getID(eEndogenous, varr), 0), recursive_variables);
                  if (d1 == Zero)
                    continue;
                  first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))] = d1;
                  tmp_derivatives.push_back(
                                            make_pair(make_pair(eq, var), make_pair(0, first_chain_rule_derivatives[make_pair(eqr, make_pair(varr, 0))])));
                }
            }
        }
      blocks_derivatives[block] = tmp_derivatives;
    }
}

void
DynamicModel::collect_block_first_order_derivatives()
{
  //! vector for an equation or a variable indicates the block number
  vector<int> equation_2_block, variable_2_block;
  unsigned int nb_blocks = getNbBlocks();
  equation_2_block = vector<int>(equation_reordered.size());
  variable_2_block = vector<int>(variable_reordered.size());
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      unsigned int block_size = getBlockSize(block);
      for (unsigned int i = 0; i < block_size; i++)
        {
          equation_2_block[getBlockEquationID(block, i)] = block;
          variable_2_block[getBlockVariableID(block, i)] = block;
        }
    }
  other_endo_block = vector<t_lag_var>(nb_blocks);
  exo_block = vector<t_lag_var>(nb_blocks);
  exo_det_block = vector<t_lag_var>(nb_blocks);
  derivative_endo = vector<t_derivative>(nb_blocks);
  derivative_other_endo = vector<t_derivative>(nb_blocks);
  derivative_exo = vector<t_derivative>(nb_blocks);
  derivative_exo_det = vector<t_derivative>(nb_blocks);
  endo_max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  other_endo_max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  exo_max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  exo_det_max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  max_leadlag_block = vector<pair<int, int> >(nb_blocks, make_pair(0, 0));
  for (first_derivatives_type::iterator it2 = first_derivatives.begin();
       it2 != first_derivatives.end(); it2++)
    {
      int eq = it2->first.first;
      int var = symbol_table.getTypeSpecificID(getSymbIDByDerivID(it2->first.second));
      int lag = getLagByDerivID(it2->first.second);
      int block_eq = equation_2_block[eq];
      int block_var = variable_2_block[var];
      t_derivative tmp_derivative;
      t_lag_var lag_var;
      switch (getTypeByDerivID(it2->first.second))
        {
        case eEndogenous:
          if (block_eq == block_var)
            {
              if (lag < 0 && lag < -endo_max_leadlag_block[block_eq].first)
                endo_max_leadlag_block[block_eq] = make_pair(-lag, endo_max_leadlag_block[block_eq].second);
              if (lag > 0 && lag > endo_max_leadlag_block[block_eq].second)
                endo_max_leadlag_block[block_eq] = make_pair(endo_max_leadlag_block[block_eq].first, lag);
              tmp_derivative = derivative_endo[block_eq];
              tmp_derivative[make_pair(lag, make_pair(eq, var))] = first_derivatives[make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, var), lag))];
              derivative_endo[block_eq] = tmp_derivative;
            }
          else
            {
              if (lag < 0 && lag < -other_endo_max_leadlag_block[block_eq].first)
                other_endo_max_leadlag_block[block_eq] = make_pair(-lag, other_endo_max_leadlag_block[block_eq].second);
              if (lag > 0 && lag > other_endo_max_leadlag_block[block_eq].second)
                other_endo_max_leadlag_block[block_eq] = make_pair(other_endo_max_leadlag_block[block_eq].first, lag);
              tmp_derivative = derivative_other_endo[block_eq];
              tmp_derivative[make_pair(lag, make_pair(eq, var))] = first_derivatives[make_pair(eq, getDerivID(symbol_table.getID(eEndogenous, var), lag))];
              derivative_other_endo[block_eq] = tmp_derivative;
              lag_var = other_endo_block[block_eq];
              if (lag_var.find(lag) == lag_var.end())
                lag_var[lag].clear();
              lag_var[lag].insert(var);
              other_endo_block[block_eq] = lag_var;
            }
          break;
        case eExogenous:
          if (lag < 0 && lag < -exo_max_leadlag_block[block_eq].first)
            exo_max_leadlag_block[block_eq] = make_pair(-lag, exo_max_leadlag_block[block_eq].second);
          if (lag > 0 && lag > exo_max_leadlag_block[block_eq].second)
            exo_max_leadlag_block[block_eq] = make_pair(exo_max_leadlag_block[block_eq].first, lag);
          tmp_derivative = derivative_exo[block_eq];
          tmp_derivative[make_pair(lag, make_pair(eq, var))] = first_derivatives[make_pair(eq, getDerivID(symbol_table.getID(eExogenous, var), lag))];
          derivative_exo[block_eq] = tmp_derivative;
          lag_var = exo_block[block_eq];
          if (lag_var.find(lag) == lag_var.end())
            lag_var[lag].clear();
          lag_var[lag].insert(var);
          exo_block[block_eq] = lag_var;
          break;
        case eExogenousDet:
          if (lag < 0 && lag < -exo_det_max_leadlag_block[block_eq].first)
            exo_det_max_leadlag_block[block_eq] = make_pair(-lag, exo_det_max_leadlag_block[block_eq].second);
          if (lag > 0 && lag > exo_det_max_leadlag_block[block_eq].second)
            exo_det_max_leadlag_block[block_eq] = make_pair(exo_det_max_leadlag_block[block_eq].first, lag);
          tmp_derivative = derivative_exo_det[block_eq];
          tmp_derivative[make_pair(lag, make_pair(eq, var))] = first_derivatives[make_pair(eq, getDerivID(symbol_table.getID(eExogenous, var), lag))];
          derivative_exo_det[block_eq] = tmp_derivative;
          lag_var = exo_det_block[block_eq];
          if (lag_var.find(lag) == lag_var.end())
            lag_var[lag].clear();
          lag_var[lag].insert(var);
          exo_det_block[block_eq] = lag_var;
          break;
        default:
          break;
        }
      if (lag < 0 && lag < -max_leadlag_block[block_eq].first)
        max_leadlag_block[block_eq] = make_pair(-lag, max_leadlag_block[block_eq].second);
      if (lag > 0 && lag > max_leadlag_block[block_eq].second)
        max_leadlag_block[block_eq] = make_pair(max_leadlag_block[block_eq].first, lag);
    }

}

void
DynamicModel::writeDynamicFile(const string &basename, bool block, bool bytecode, bool use_dll, int order) const
{
  int r;
  if (block && bytecode)
    writeModelEquationsCodeOrdered(basename + "_dynamic", basename, map_idx);
  else if (block && !bytecode)
    {
#ifdef _WIN32
      r = mkdir(basename.c_str());
#else
      r = mkdir(basename.c_str(), 0777);
#endif
      if (r < 0 && errno != EEXIST)
        {
          perror("ERROR");
          exit(EXIT_FAILURE);
        }
      writeSparseDynamicMFile(basename + "_dynamic", basename);
    }
  else if (use_dll)
    writeDynamicCFile(basename + "_dynamic", order);
  else
    writeDynamicMFile(basename + "_dynamic");
}

void
DynamicModel::toStatic(StaticModel &static_model) const
{
  // Convert model local variables (need to be done first)
  for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    static_model.AddLocalVariable(symbol_table.getName(it->first), it->second->toStatic(static_model));

  // Convert equations
  for (vector<BinaryOpNode *>::const_iterator it = equations.begin();
       it != equations.end(); it++)
    static_model.addEquation((*it)->toStatic(static_model));

  // Convert auxiliary equations
  for (deque<BinaryOpNode *>::const_iterator it = aux_equations.begin();
       it != aux_equations.end(); it++)
    static_model.addAuxEquation((*it)->toStatic(static_model));
}

void
DynamicModel::computeDerivIDs()
{
  set<pair<int, int> > dynvars;

  for (int i = 0; i < (int) equations.size(); i++)
    equations[i]->collectVariables(eEndogenous, dynvars);

  dynJacobianColsNbr = dynvars.size();

  for (int i = 0; i < (int) equations.size(); i++)
    {
      equations[i]->collectVariables(eExogenous, dynvars);
      equations[i]->collectVariables(eExogenousDet, dynvars);
      equations[i]->collectVariables(eParameter, dynvars);
    }

  for (set<pair<int, int> >::const_iterator it = dynvars.begin();
       it != dynvars.end(); it++)
    {
      int lag = it->second;
      SymbolType type = symbol_table.getType(it->first);

      /* Setting maximum and minimum lags.

         We don't want these to be affected by lead/lags on parameters: they
         are accepted for facilitating variable flipping, but are simply
         ignored. */
      if (max_lead < lag && type != eParameter)
        max_lead = lag;
      else if (-max_lag > lag && type != eParameter)
        max_lag = -lag;

      switch (type)
        {
        case eEndogenous:
          if (max_endo_lead < lag)
            max_endo_lead = lag;
          else if (-max_endo_lag > lag)
            max_endo_lag = -lag;
          break;
        case eExogenous:
          if (max_exo_lead < lag)
            max_exo_lead = lag;
          else if (-max_exo_lag > lag)
            max_exo_lag = -lag;
          break;
        case eExogenousDet:
          if (max_exo_det_lead < lag)
            max_exo_det_lead = lag;
          else if (-max_exo_det_lag > lag)
            max_exo_det_lag = -lag;
          break;
        default:
          break;
        }

      // Create a new deriv_id
      int deriv_id = deriv_id_table.size();

      deriv_id_table[*it] = deriv_id;
      inv_deriv_id_table.push_back(*it);
    }
}

SymbolType
DynamicModel::getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  return symbol_table.getType(getSymbIDByDerivID(deriv_id));
}

int
DynamicModel::getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].second;
}

int
DynamicModel::getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  if (deriv_id < 0 || deriv_id >= (int) inv_deriv_id_table.size())
    throw UnknownDerivIDException();

  return inv_deriv_id_table[deriv_id].first;
}

int
DynamicModel::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  deriv_id_table_t::const_iterator it = deriv_id_table.find(make_pair(symb_id, lag));
  if (it == deriv_id_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::computeDynJacobianCols(bool jacobianExo)
{
  /* Sort the dynamic endogenous variables by lexicographic order over (lag, type_specific_symbol_id)
     and fill the dynamic columns for exogenous and exogenous deterministic */
  map<pair<int, int>, int> ordered_dyn_endo;

  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      const int &symb_id = it->first.first;
      const int &lag = it->first.second;
      const int &deriv_id = it->second;
      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id);

      switch (type)
        {
        case eEndogenous:
          ordered_dyn_endo[make_pair(lag, tsid)] = deriv_id;
          break;
        case eExogenous:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + tsid;
          break;
        case eExogenousDet:
          // At this point, dynJacobianColsNbr contains the number of dynamic endogenous
          if (jacobianExo)
            dyn_jacobian_cols_table[deriv_id] = dynJacobianColsNbr + symbol_table.exo_nbr() + tsid;
          break;
        case eParameter:
          // We don't assign a dynamic jacobian column to parameters
          break;
        default:
          // Shut up GCC
          cerr << "DynamicModel::computeDynJacobianCols: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
    }

  // Fill in dynamic jacobian columns for endogenous
  int sorted_id = 0;
  for (map<pair<int, int>, int>::const_iterator it = ordered_dyn_endo.begin();
       it != ordered_dyn_endo.end(); it++)
    dyn_jacobian_cols_table[it->second] = sorted_id++;

  // Set final value for dynJacobianColsNbr
  if (jacobianExo)
    dynJacobianColsNbr += symbol_table.exo_nbr() + symbol_table.exo_det_nbr();
}

int
DynamicModel::getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException)
{
  map<int, int>::const_iterator it = dyn_jacobian_cols_table.find(deriv_id);
  if (it == dyn_jacobian_cols_table.end())
    throw UnknownDerivIDException();
  else
    return it->second;
}

void
DynamicModel::computeParamsDerivatives()
{
  for (deriv_id_table_t::const_iterator it = deriv_id_table.begin();
       it != deriv_id_table.end(); it++)
    {
      if (symbol_table.getType(it->first.first) != eParameter)
        continue;

      int param = it->second;

      for (int eq = 0; eq < (int) equations.size(); eq++)
        {
          NodeID d1 = equations[eq]->getDerivative(param);
          if (d1 == Zero)
            continue;
          residuals_params_derivatives[make_pair(eq, param)] = d1;
        }

      for (first_derivatives_type::const_iterator it2 = first_derivatives.begin();
           it2 != first_derivatives.end(); it2++)
        {
          int eq = it2->first.first;
          int var = it2->first.second;
          NodeID d1 = it2->second;

          NodeID d2 = d1->getDerivative(param);
          if (d2 == Zero)
            continue;
          jacobian_params_derivatives[make_pair(eq, make_pair(var, param))] = d2;
        }
    }
}

void
DynamicModel::computeParamsDerivativesTemporaryTerms()
{
  map<NodeID, int> reference_count;
  params_derivs_temporary_terms.clear();

  for (first_derivatives_type::iterator it = residuals_params_derivatives.begin();
       it != residuals_params_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, params_derivs_temporary_terms, true);

  for (second_derivatives_type::iterator it = jacobian_params_derivatives.begin();
       it != jacobian_params_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, params_derivs_temporary_terms, true);
}

void
DynamicModel::writeParamsDerivativesFile(const string &basename) const
{
  if (!residuals_params_derivatives.size()
      && !jacobian_params_derivatives.size())
    return;

  string filename = basename + "_params_derivs.m";

  ofstream paramsDerivsFile;
  paramsDerivsFile.open(filename.c_str(), ios::out | ios::binary);
  if (!paramsDerivsFile.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }
  paramsDerivsFile << "function [rp, gp] = " << basename << "_params_derivs(y, x, params, it_)" << endl
                   << "%" << endl
                   << "% Warning : this file is generated automatically by Dynare" << endl
                   << "%           from model file (.mod)" << endl << endl;

  if (isUnaryOpUsed(oSteadyState))
    paramsDerivsFile << "global oo_;" << endl << endl;

  writeTemporaryTerms(params_derivs_temporary_terms, paramsDerivsFile, oMatlabDynamicModel);

  // Write parameter derivative
  paramsDerivsFile << "rp = zeros(" << equation_number() << ", "
                   << symbol_table.param_nbr() << ");" << endl;

  for (first_derivatives_type::const_iterator it = residuals_params_derivatives.begin();
       it != residuals_params_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int param = it->first.second;
      NodeID d1 = it->second;

      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      paramsDerivsFile << "rp(" << eq+1 << ", " << param_col << ") = ";
      d1->writeOutput(paramsDerivsFile, oMatlabDynamicModel, params_derivs_temporary_terms);
      paramsDerivsFile << ";" << endl;
    }

  // Write jacobian derivatives
  paramsDerivsFile << "gp = zeros(" << equation_number() << ", " << dynJacobianColsNbr << ", "
                   << symbol_table.param_nbr() << ");" << endl;

  for (second_derivatives_type::const_iterator it = jacobian_params_derivatives.begin();
       it != jacobian_params_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var = it->first.second.first;
      int param = it->first.second.second;
      NodeID d2 = it->second;

      int var_col = getDynJacobianCol(var) + 1;
      int param_col = symbol_table.getTypeSpecificID(getSymbIDByDerivID(param)) + 1;

      paramsDerivsFile << "gp(" << eq+1 << ", " << var_col << ", " << param_col << ") = ";
      d2->writeOutput(paramsDerivsFile, oMatlabDynamicModel, params_derivs_temporary_terms);
      paramsDerivsFile << ";" << endl;
    }

  paramsDerivsFile.close();
}

void
DynamicModel::writeChainRuleDerivative(ostream &output, int eqr, int varr, int lag,
                                       ExprNodeOutputType output_type,
                                       const temporary_terms_type &temporary_terms) const
{
  map<pair<int, pair<int, int> >, NodeID>::const_iterator it = first_chain_rule_derivatives.find(make_pair(eqr, make_pair(varr, lag)));
  if (it != first_chain_rule_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
DynamicModel::writeLatexFile(const string &basename) const
{
  writeLatexModelFile(basename + "_dynamic.tex", oLatexDynamicModel);
}

void
DynamicModel::jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << eq_nb + 1 << "," << col_nb + 1;
  else
    output << eq_nb + col_nb *equations.size();
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
DynamicModel::sparseHelper(int order, ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const
{
  output << "v" << order << LEFT_ARRAY_SUBSCRIPT(output_type);
  if (IS_MATLAB(output_type))
    output << row_nb + 1 << "," << col_nb + 1;
  else
    output << row_nb + col_nb * NNZDerivatives[order-1];
  output << RIGHT_ARRAY_SUBSCRIPT(output_type);
}

void
DynamicModel::substituteEndoLeadGreaterThanTwo()
{
  substituteLeadLagInternal(avEndoLead);
}

void
DynamicModel::substituteEndoLagGreaterThanTwo()
{
  substituteLeadLagInternal(avEndoLag);
}

void
DynamicModel::substituteExoLead()
{
  substituteLeadLagInternal(avExoLead);
}

void
DynamicModel::substituteExoLag()
{
  substituteLeadLagInternal(avExoLag);
}

void
DynamicModel::substituteLeadLagInternal(aux_var_t type)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in used model local variables
  set<int> used_local_vars;
  for (size_t i = 0; i < equations.size(); i++)
    equations[i]->collectModelLocalVariables(used_local_vars);

  for (set<int>::const_iterator it = used_local_vars.begin();
       it != used_local_vars.end(); ++it)
    {
      const NodeID value = local_variables_table.find(*it)->second;
      NodeID subst;
      switch (type)
        {
        case avEndoLead:
          subst = value->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          break;
        case avEndoLag:
          subst = value->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case avExoLead:
          subst = value->substituteExoLead(subst_table, neweqs);
          break;
        case avExoLag:
          subst = value->substituteExoLag(subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      local_variables_table[*it] = subst;
    }

  // Substitute in equations
  for (int i = 0; i < (int) equations.size(); i++)
    {
      NodeID subst;
      switch (type)
        {
        case avEndoLead:
          subst = equations[i]->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          break;
        case avEndoLag:
          subst = equations[i]->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
          break;
        case avExoLead:
          subst = equations[i]->substituteExoLead(subst_table, neweqs);
          break;
        case avExoLag:
          subst = equations[i]->substituteExoLag(subst_table, neweqs);
          break;
        default:
          cerr << "DynamicModel::substituteLeadLagInternal: impossible case" << endl;
          exit(EXIT_FAILURE);
        }
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(subst);
      assert(substeq != NULL);
      equations[i] = substeq;
    }

  // Add new equations
  for (int i = 0; i < (int) neweqs.size(); i++)
    addEquation(neweqs[i]);

  // Add the new set of equations at the *beginning* of aux_equations
  copy(neweqs.rbegin(), neweqs.rend(), front_inserter(aux_equations));

  if (neweqs.size() > 0)
    {
      cout << "Substitution of ";
      switch (type)
        {
        case avEndoLead:
          cout << "endo leads >= 2";
          break;
        case avEndoLag:
          cout << "endo lags >= 2";
          break;
        case avExoLead:
          cout << "exo leads";
          break;
        case avExoLag:
          cout << "exo lags";
          break;
        case avExpectation:
          cout << "expectation";
          break;
        }
      cout << ": added " << neweqs.size() << " auxiliary variables and equations." << endl;
    }
}

void
DynamicModel::substituteExpectation(bool partial_information_model)
{
  ExprNode::subst_table_t subst_table;
  vector<BinaryOpNode *> neweqs;

  // Substitute in model local variables
  for (map<int, NodeID>::iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    it->second = it->second->substituteExpectation(subst_table, neweqs, partial_information_model);

  // Substitute in equations
  for (int i = 0; i < (int) equations.size(); i++)
    {
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(equations[i]->substituteExpectation(subst_table, neweqs, partial_information_model));
      assert(substeq != NULL);
      equations[i] = substeq;
    }

  // Add new equations
  for (int i = 0; i < (int) neweqs.size(); i++)
    addEquation(neweqs[i]);

  // Add the new set of equations at the *beginning* of aux_equations
  copy(neweqs.rbegin(), neweqs.rend(), front_inserter(aux_equations));

  if (neweqs.size() > 0)
    {
      if (partial_information_model)
        cout << "Substitution of Expectation operator: added " << subst_table.size() << " auxiliary variables and " << neweqs.size() << " auxiliary equations." << endl;
      else
        cout << "Substitution of Expectation operator: added " << neweqs.size() << " auxiliary variables and equations." << endl;
    }
}

void
DynamicModel::transformPredeterminedVariables()
{
  for (int i = 0; i < (int) equations.size(); i++)
    {
      BinaryOpNode *substeq = dynamic_cast<BinaryOpNode *>(equations[i]->decreaseLeadsLagsPredeterminedVariables());
      assert(substeq != NULL);
      equations[i] = substeq;
    }
}

void
DynamicModel::fillEvalContext(eval_context_type &eval_context) const
{
  // First, auxiliary variables
  for (deque<BinaryOpNode *>::const_iterator it = aux_equations.begin();
       it != aux_equations.end(); it++)
    {
      assert((*it)->get_op_code() == oEqual);
      VariableNode *auxvar = dynamic_cast<VariableNode *>((*it)->get_arg1());
      assert(auxvar != NULL);
      try
        {
          double val = (*it)->get_arg2()->eval(eval_context);
          eval_context[auxvar->get_symb_id()] = val;
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }

  // Second, model local variables
  for (map<int, NodeID>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    {
      try
        {
          const NodeID expression = it->second;
          double val = expression->eval(eval_context);
          eval_context[it->first] = val;
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }
}

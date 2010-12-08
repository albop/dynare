/*
 * Copyright (C) 2007-2009 Dynare Team
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
#include <iterator>
#include <algorithm>

// For select1st()
#ifdef __GNUC__
# include <ext/functional>
using namespace __gnu_cxx;
#endif

#include <cassert>
#include <cmath>

#include "ExprNode.hh"
#include "DataTree.hh"
#include "ModFile.hh"

ExprNode::ExprNode(DataTree &datatree_arg) : datatree(datatree_arg), preparedForDerivation(false)
{
  // Add myself to datatree
  datatree.node_list.push_back(this);

  // Set my index and increment counter
  idx = datatree.node_counter++;
}

ExprNode::~ExprNode()
{
}

NodeID
ExprNode::getDerivative(int deriv_id)
{
  if (!preparedForDerivation)
    prepareForDerivation();

  // Return zero if derivative is necessarily null (using symbolic a priori)
  set<int>::const_iterator it = non_null_derivatives.find(deriv_id);
  if (it == non_null_derivatives.end())
    return datatree.Zero;

  // If derivative is stored in cache, use the cached value, otherwise compute it (and cache it)
  map<int, NodeID>::const_iterator it2 = derivatives.find(deriv_id);
  if (it2 != derivatives.end())
    return it2->second;
  else
    {
      NodeID d = computeDerivative(deriv_id);
      derivatives[deriv_id] = d;
      return d;
    }
}

int
ExprNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
{
  // For a constant, a variable, or a unary op, the precedence is maximal
  return 100;
}

int
ExprNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

void
ExprNode::collectEndogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eEndogenous, symb_ids);
  for (set<pair<int, int> >::const_iterator it = symb_ids.begin();
       it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::collectExogenous(set<pair<int, int> > &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eExogenous, symb_ids);
  for (set<pair<int, int> >::const_iterator it = symb_ids.begin();
       it != symb_ids.end(); it++)
    result.insert(make_pair(datatree.symbol_table.getTypeSpecificID(it->first), it->second));
}

void
ExprNode::collectModelLocalVariables(set<int> &result) const
{
  set<pair<int, int> > symb_ids;
  collectVariables(eModelLocalVariable, symb_ids);
  transform(symb_ids.begin(), symb_ids.end(), inserter(result, result.begin()),
            select1st<pair<int, int> >());
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                bool is_matlab) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                map<NodeID, pair<int, int> > &first_occurence,
                                int Curr_block,
                                vector<vector<temporary_terms_type> > &v_temporary_terms,
                                int equation) const
{
  // Nothing to do for a terminal node
}

pair<int, NodeID >
ExprNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  return (make_pair(0, (NodeID) NULL));
}

void
ExprNode::writeOutput(ostream &output)
{
  writeOutput(output, oMatlabOutsideModel, temporary_terms_type());
}

VariableNode *
ExprNode::createEndoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  int n = maxEndoLead();
  assert(n >= 2);

  subst_table_t::const_iterator it = subst_table.find(this);
  if (it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  NodeID substexpr = decreaseLeadsLags(n-1);
  int lag = n-2;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      NodeID orig_expr = decreaseLeadsLags(lag);
      it = subst_table.find(orig_expr);
      if (it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addEndoLeadAuxiliaryVar(orig_expr->idx);
          neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr)));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr) != NULL);
          subst_table[orig_expr] = dynamic_cast<VariableNode *>(substexpr);
        }
      else
        substexpr = const_cast<VariableNode *>(it->second);

      lag--;
    }

  return dynamic_cast<VariableNode *>(substexpr);
}

VariableNode *
ExprNode::createExoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  int n = maxExoLead();
  assert(n >= 1);

  subst_table_t::const_iterator it = subst_table.find(this);
  if (it != subst_table.end())
    return const_cast<VariableNode *>(it->second);

  NodeID substexpr = decreaseLeadsLags(n);
  int lag = n-1;

  // Each iteration tries to create an auxvar such that auxvar(+1)=expr(-lag)
  // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to expr(-lag-1) (resp. expr(-lag))
  while (lag >= 0)
    {
      NodeID orig_expr = decreaseLeadsLags(lag);
      it = subst_table.find(orig_expr);
      if (it == subst_table.end())
        {
          int symb_id = datatree.symbol_table.addExoLeadAuxiliaryVar(orig_expr->idx);
          neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(symb_id, 0), substexpr)));
          substexpr = datatree.AddVariable(symb_id, +1);
          assert(dynamic_cast<VariableNode *>(substexpr) != NULL);
          subst_table[orig_expr] = dynamic_cast<VariableNode *>(substexpr);
        }
      else
        substexpr = const_cast<VariableNode *>(it->second);

      lag--;
    }

  return dynamic_cast<VariableNode *>(substexpr);
}

bool
ExprNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
ExprNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

NumConstNode::NumConstNode(DataTree &datatree_arg, int id_arg) :
  ExprNode(datatree_arg),
  id(id_arg)
{
  // Add myself to the num const map
  datatree.num_const_node_map[id] = this;
}

void
NumConstNode::prepareForDerivation()
{
  preparedForDerivation = true;
  // All derivatives are null, so non_null_derivatives is left empty
}

NodeID
NumConstNode::computeDerivative(int deriv_id)
{
  return datatree.Zero;
}

void
NumConstNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
}

void
NumConstNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    if (output_type == oMatlabDynamicModelSparse)
      output << "T" << idx << "(it_)";
    else
      output << "T" << idx;
  else
    output << datatree.num_constants.get(id);
}

double
NumConstNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  return (datatree.num_constants.getDouble(id));
}

void
NumConstNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  FLDC_ fldc(datatree.num_constants.getDouble(id));
  fldc.write(CompileCode);
}

void
NumConstNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
}

pair<int, NodeID >
NumConstNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  return (make_pair(0, datatree.AddNumConstant(datatree.num_constants.get(id))));
}

NodeID
NumConstNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  return datatree.Zero;
}

NodeID
NumConstNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddNumConstant(datatree.num_constants.get(id));
}

int
NumConstNode::maxEndoLead() const
{
  return 0;
}

int
NumConstNode::maxExoLead() const
{
  return 0;
}

int
NumConstNode::maxEndoLag() const
{
  return 0;
}

int
NumConstNode::maxExoLag() const
{
  return 0;
}

NodeID
NumConstNode::decreaseLeadsLags(int n) const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::decreaseLeadsLagsPredeterminedVariables() const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  return const_cast<NumConstNode *>(this);
}

NodeID
NumConstNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<NumConstNode *>(this);
}

bool
NumConstNode::isNumConstNodeEqualTo(double value) const
{
  if (datatree.num_constants.getDouble(id) == value)
    return true;
  else
    return false;
}

bool
NumConstNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

VariableNode::VariableNode(DataTree &datatree_arg, int symb_id_arg, int lag_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  type(datatree.symbol_table.getType(symb_id_arg)),
  lag(lag_arg)
{
  // Add myself to the variable map
  datatree.variable_node_map[make_pair(symb_id, lag)] = this;

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  assert(type != eUnknownFunction
         && (lag == 0 || (type != eModelLocalVariable && type != eModFileLocalVariable)));
}

void
VariableNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  // Fill in non_null_derivatives
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      // For a variable or a parameter, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(datatree.getDerivID(symb_id, lag));
      break;
    case eModelLocalVariable:
      datatree.local_variables_table[symb_id]->prepareForDerivation();
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.local_variables_table[symb_id]->non_null_derivatives;
      break;
    case eModFileLocalVariable:
      // Such a variable is never derived
      break;
    case eUnknownFunction:
      cerr << "VariableNode::prepareForDerivation: impossible case" << endl;
      exit(EXIT_FAILURE);
    }
}

NodeID
VariableNode::computeDerivative(int deriv_id)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      else
        return datatree.Zero;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getDerivative(deriv_id);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnknownFunction:
      cerr << "Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
VariableNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
}

void
VariableNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  if (IS_LATEX(output_type))
    {
      if (output_type == oLatexDynamicSteadyStateOperator)
        output << "\\bar{";
      output << datatree.symbol_table.getTeXName(symb_id);
      if (output_type == oLatexDynamicModel
          && (type == eEndogenous || type == eExogenous || type == eExogenousDet || type == eModelLocalVariable))
        {
          output << "_{t";
          if (lag != 0)
            {
              if (lag > 0)
                output << "+";
              output << lag;
            }
          output << "}";
        }
      else if (output_type == oLatexDynamicSteadyStateOperator)
        output << "}";
      return;
    }

  int i;
  int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
  switch (type)
    {
    case eParameter:
      if (output_type == oMatlabOutsideModel)
        output << "M_.params" << "(" << tsid + 1 << ")";
      else
        output << "params" << LEFT_ARRAY_SUBSCRIPT(output_type) << tsid + ARRAY_SUBSCRIPT_OFFSET(output_type) << RIGHT_ARRAY_SUBSCRIPT(output_type);
      break;

    case eModelLocalVariable:
    case eModFileLocalVariable:
      if (output_type == oMatlabDynamicModelSparse || output_type == oMatlabStaticModelSparse || output_type == oMatlabDynamicModelSparseLocalTemporaryTerms)
        {
          output << "(";
          datatree.local_variables_table[symb_id]->writeOutput(output, output_type, temporary_terms);
          output << ")";
        }
      else
        output << datatree.symbol_table.getName(symb_id);
      break;

    case eEndogenous:
      switch (output_type)
        {
        case oMatlabDynamicModel:
        case oCDynamicModel:
          i = datatree.getDynJacobianCol(datatree.getDerivID(symb_id, lag)) + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
          output <<  "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabDynamicModelSparse:
        case oMatlabDynamicModelSparseLocalTemporaryTerms:
          i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
          if (lag > 0)
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_+" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else if (lag < 0)
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_" << lag << ", " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          else
            output << "y" << LEFT_ARRAY_SUBSCRIPT(output_type) << "it_, " << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          output << "oo_.steady_state(" << tsid + 1 << ")";
          break;
        case oMatlabDynamicSteadyStateOperator:
        case oMatlabDynamicSparseSteadyStateOperator:
          output << "oo_.steady_state(" << tsid + 1 << ")";
          break;
        case oCDynamicSteadyStateOperator:
          output << "steady_state[" << tsid << "]";
          break;
        default:
          assert(false);
        }
      break;

    case eExogenous:
      i = tsid + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case oMatlabDynamicModel:
        case oMatlabDynamicModelSparse:
        case oMatlabDynamicModelSparseLocalTemporaryTerms:
          if (lag > 0)
            output <<  "x(it_+" << lag << ", " << i << ")";
          else if (lag < 0)
            output <<  "x(it_" << lag << ", " << i << ")";
          else
            output <<  "x(it_, " << i << ")";
          break;
        case oCDynamicModel:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        case oMatlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_steady_state(" << i << ")";
          break;
        default:
          assert(false);
        }
      break;

    case eExogenousDet:
      i = tsid + datatree.symbol_table.exo_nbr() + ARRAY_SUBSCRIPT_OFFSET(output_type);
      switch (output_type)
        {
        case oMatlabDynamicModel:
        case oMatlabDynamicModelSparse:
        case oMatlabDynamicModelSparseLocalTemporaryTerms:
          if (lag > 0)
            output <<  "x(it_+" << lag << ", " << i << ")";
          else if (lag < 0)
            output <<  "x(it_" << lag << ", " << i << ")";
          else
            output <<  "x(it_, " << i << ")";
          break;
        case oCDynamicModel:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case oMatlabStaticModel:
        case oMatlabStaticModelSparse:
          output << "x" << LEFT_ARRAY_SUBSCRIPT(output_type) << i << RIGHT_ARRAY_SUBSCRIPT(output_type);
          break;
        case oMatlabOutsideModel:
          assert(lag == 0);
          output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
          break;
        case oMatlabDynamicSteadyStateOperator:
          output <<  "oo_.exo_det_steady_state(" << tsid + 1 << ")";
          break;
        default:
          assert(false);
        }
      break;

    case eUnknownFunction:
      cerr << "Impossible case" << endl;
      exit(EXIT_FAILURE);
    }
}

double
VariableNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  eval_context_type::const_iterator it = eval_context.find(symb_id);
  if (it == eval_context.end())
    throw EvalException();

  return it->second;
}

void
VariableNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  if (type == eModelLocalVariable || type == eModFileLocalVariable)
    datatree.local_variables_table[symb_id]->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  else
    {
      int tsid = datatree.symbol_table.getTypeSpecificID(symb_id);
      if (type == eExogenousDet)
        tsid += datatree.symbol_table.exo_nbr();
      if (!lhs_rhs)
        {
          if (dynamic)
            {
              if (steady_dynamic)  // steady state values in a dynamic model
                {
                  FLDVS_ fldvs(type, tsid);
                  fldvs.write(CompileCode);
                }
              else
                {
                  if (type == eParameter)
                    {
                      FLDV_ fldv(type, tsid);
                      fldv.write(CompileCode);
                    }
                  else
                    {
                      FLDV_ fldv(type, tsid, lag);
                      fldv.write(CompileCode);
                    }
                }
            }
          else
            {
              FLDSV_ fldsv(type, tsid);
              fldsv.write(CompileCode);
            }
        }
      else
        {
          if (dynamic)
            {
              if (steady_dynamic)  // steady state values in a dynamic model
                {
                  cerr << "Impossible case: steady_state in rhs of equation" << endl;
                  exit(EXIT_FAILURE);
                }
              else
                {
                  if (type == eParameter)
                    {
                      FSTPV_ fstpv(type, tsid);
                      fstpv.write(CompileCode);
                    }
                  else
                    {
                      FSTPV_ fstpv(type, tsid, lag);
                      fstpv.write(CompileCode);
                    }
                }
            }
          else
            {
              FSTPSV_ fstpsv(type, tsid);
              fstpsv.write(CompileCode);
            }
        }
    }
}

void
VariableNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    map<NodeID, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    vector<vector<temporary_terms_type> > &v_temporary_terms,
                                    int equation) const
{
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
}

void
VariableNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  if (type == type_arg)
    result.insert(make_pair(symb_id, lag));
  if (type == eModelLocalVariable)
    datatree.local_variables_table[symb_id]->collectVariables(type_arg, result);
}

pair<int, NodeID>
VariableNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  if (type == eEndogenous)
    {
      if (datatree.symbol_table.getTypeSpecificID(symb_id) == var_endo && lag == 0)
        return (make_pair(1, (NodeID) NULL));
      else
        return (make_pair(0, datatree.AddVariableInternal(symb_id, lag)));
    }
  else
    {
      if (type == eParameter)
        return (make_pair(0, datatree.AddVariableInternal(symb_id, 0)));
      else
        return (make_pair(0, datatree.AddVariableInternal(symb_id, lag)));
    }
}

NodeID
VariableNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eParameter:
      if (deriv_id == datatree.getDerivID(symb_id, lag))
        return datatree.One;
      else
        {
          //if there is in the equation a recursive variable we could use a chaine rule derivation
          map<int, NodeID>::const_iterator it = recursive_variables.find(datatree.getDerivID(symb_id, lag));
          if (it != recursive_variables.end())
            {
              map<int, NodeID>::const_iterator it2 = derivatives.find(deriv_id);
              if (it2 != derivatives.end())
                return it2->second;
              else
                {
                  map<int, NodeID> recursive_vars2(recursive_variables);
                  recursive_vars2.erase(it->first);
                  //NodeID c = datatree.AddNumConstant("1");
                  NodeID d = datatree.AddUMinus(it->second->getChainRuleDerivative(deriv_id, recursive_vars2));
                  //d = datatree.AddTimes(c, d);
                  derivatives[deriv_id] = d;
                  return d;
                }
            }
          else
            return datatree.Zero;
        }
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getChainRuleDerivative(deriv_id, recursive_variables);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(EXIT_FAILURE);
    case eUnknownFunction:
      cerr << "Impossible case!" << endl;
      exit(EXIT_FAILURE);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
VariableNode::toStatic(DataTree &static_datatree) const
{
  return static_datatree.AddVariable(datatree.symbol_table.getName(symb_id));
}

int
VariableNode::maxEndoLead() const
{
  switch (type)
    {
    case eEndogenous:
      return max(lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxEndoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLead() const
{
  switch (type)
    {
    case eExogenous:
      return max(lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxExoLead();
    default:
      return 0;
    }
}

int
VariableNode::maxEndoLag() const
{
  switch (type)
    {
    case eEndogenous:
      return max(-lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxEndoLag();
    default:
      return 0;
    }
}

int
VariableNode::maxExoLag() const
{
  switch (type)
    {
    case eExogenous:
      return max(-lag, 0);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->maxExoLag();
    default:
      return 0;
    }
}

NodeID
VariableNode::decreaseLeadsLags(int n) const
{
  switch (type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
      return datatree.AddVariable(symb_id, lag-n);
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->decreaseLeadsLags(n);
    default:
      return const_cast<VariableNode *>(this);
    }
}

NodeID
VariableNode::decreaseLeadsLagsPredeterminedVariables() const
{
  if (datatree.symbol_table.isPredetermined(symb_id))
    return decreaseLeadsLags(1);
  else
    return const_cast<VariableNode *>(this);
}

NodeID
VariableNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID value;
  switch (type)
    {
    case eEndogenous:
      if (lag <= 1)
        return const_cast<VariableNode *>(this);
      else
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxEndoLead() <= 1)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

NodeID
VariableNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  VariableNode *substexpr;
  NodeID value;
  subst_table_t::const_iterator it;
  int cur_lag;
  switch (type)
    {
    case eEndogenous:
      if (lag >= -1)
        return const_cast<VariableNode *>(this);

      it = subst_table.find(this);
      if (it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, -1);
      cur_lag = -2;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          it = subst_table.find(orig_expr);
          if (it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addEndoLagAuxiliaryVar(symb_id, cur_lag+1);
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr)));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxEndoLag() <= 1)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

NodeID
VariableNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID value;
  switch (type)
    {
    case eExogenous:
      if (lag <= 0)
        return const_cast<VariableNode *>(this);
      else
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxExoLead() == 0)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteExoLead(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

NodeID
VariableNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  VariableNode *substexpr;
  NodeID value;
  subst_table_t::const_iterator it;
  int cur_lag;
  switch (type)
    {
    case eExogenous:
      if (lag >= 0)
        return const_cast<VariableNode *>(this);

      it = subst_table.find(this);
      if (it != subst_table.end())
        return const_cast<VariableNode *>(it->second);

      substexpr = datatree.AddVariable(symb_id, 0);
      cur_lag = -1;

      // Each iteration tries to create an auxvar such that auxvar(-1)=curvar(cur_lag)
      // At the beginning (resp. end) of each iteration, substexpr is an expression (possibly an auxvar) equivalent to curvar(cur_lag+1) (resp. curvar(cur_lag))
      while (cur_lag >= lag)
        {
          VariableNode *orig_expr = datatree.AddVariable(symb_id, cur_lag);
          it = subst_table.find(orig_expr);
          if (it == subst_table.end())
            {
              int aux_symb_id = datatree.symbol_table.addExoLagAuxiliaryVar(symb_id, cur_lag+1);
              neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(datatree.AddVariable(aux_symb_id, 0), substexpr)));
              substexpr = datatree.AddVariable(aux_symb_id, -1);
              subst_table[orig_expr] = substexpr;
            }
          else
            substexpr = const_cast<VariableNode *>(it->second);

          cur_lag--;
        }
      return substexpr;

    case eModelLocalVariable:
      value = datatree.local_variables_table[symb_id];
      if (value->maxExoLag() == 0)
        return const_cast<VariableNode *>(this);
      else
        return value->substituteExoLag(subst_table, neweqs);
    default:
      return const_cast<VariableNode *>(this);
    }
}

NodeID
VariableNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  return const_cast<VariableNode *>(this);
}

bool
VariableNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
VariableNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  if (type == type_arg && datatree.symbol_table.getTypeSpecificID(symb_id) == variable_id && lag == lag_arg)
    return true;
  else
    return false;
}

UnaryOpNode::UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg, const int expectation_information_set_arg) :
  ExprNode(datatree_arg),
  arg(arg_arg),
  expectation_information_set(expectation_information_set_arg),
  op_code(op_code_arg)
{
  // Add myself to the unary op map
  datatree.unary_op_node_map[make_pair(make_pair(arg, op_code),
                                       expectation_information_set)] = this;
}

void
UnaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg->prepareForDerivation();

  // Non-null derivatives are those of the argument
  non_null_derivatives = arg->non_null_derivatives;
}

NodeID
UnaryOpNode::composeDerivatives(NodeID darg)
{
  NodeID t11, t12, t13;

  switch (op_code)
    {
    case oUminus:
      return datatree.AddUMinus(darg);
    case oExp:
      return datatree.AddTimes(darg, this);
    case oLog:
      return datatree.AddDivide(darg, arg);
    case oLog10:
      t11 = datatree.AddExp(datatree.One);
      t12 = datatree.AddLog10(t11);
      t13 = datatree.AddDivide(darg, arg);
      return datatree.AddTimes(t12, t13);
    case oCos:
      t11 = datatree.AddSin(arg);
      t12 = datatree.AddUMinus(t11);
      return datatree.AddTimes(darg, t12);
    case oSin:
      t11 = datatree.AddCos(arg);
      return datatree.AddTimes(darg, t11);
    case oTan:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddPlus(t11, datatree.One);
      return datatree.AddTimes(darg, t12);
    case oAcos:
      t11 = datatree.AddSin(this);
      t12 = datatree.AddDivide(darg, t11);
      return datatree.AddUMinus(t12);
    case oAsin:
      t11 = datatree.AddCos(this);
      return datatree.AddDivide(darg, t11);
    case oAtan:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddPlus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case oCosh:
      t11 = datatree.AddSinh(arg);
      return datatree.AddTimes(darg, t11);
    case oSinh:
      t11 = datatree.AddCosh(arg);
      return datatree.AddTimes(darg, t11);
    case oTanh:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oAcosh:
      t11 = datatree.AddSinh(this);
      return datatree.AddDivide(darg, t11);
    case oAsinh:
      t11 = datatree.AddCosh(this);
      return datatree.AddDivide(darg, t11);
    case oAtanh:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oSqrt:
      t11 = datatree.AddPlus(this, this);
      return datatree.AddDivide(darg, t11);
    case oSteadyState:
      if (datatree.isDynamic())
        return datatree.Zero;
      else
        return darg;
    case oExpectation:
      assert(0);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
UnaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg = arg->getDerivative(deriv_id);
  return composeDerivatives(darg);
}

int
UnaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    return 0;

  int cost = arg->cost(temporary_terms, is_matlab);

  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oUminus:
        return cost + 70;
      case oExp:
        return cost + 160;
      case oLog:
        return cost + 300;
      case oLog10:
        return cost + 16000;
      case oCos:
      case oSin:
      case oCosh:
        return cost + 210;
      case oTan:
        return cost + 230;
      case oAcos:
        return cost + 300;
      case oAsin:
        return cost + 310;
      case oAtan:
        return cost + 140;
      case oSinh:
        return cost + 240;
      case oTanh:
        return cost + 190;
      case oAcosh:
        return cost + 770;
      case oAsinh:
        return cost + 460;
      case oAtanh:
        return cost + 350;
      case oSqrt:
        return cost + 570;
      case oSteadyState:
      case oExpectation:
        return cost;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oUminus:
        return cost + 3;
      case oExp:
      case oAcosh:
        return cost + 210;
      case oLog:
        return cost + 137;
      case oLog10:
        return cost + 139;
      case oCos:
      case oSin:
        return cost + 160;
      case oTan:
        return cost + 170;
      case oAcos:
      case oAtan:
        return cost + 190;
      case oAsin:
        return cost + 180;
      case oCosh:
      case oSinh:
      case oTanh:
        return cost + 240;
      case oAsinh:
        return cost + 220;
      case oAtanh:
        return cost + 150;
      case oSqrt:
        return cost + 90;
      case oSteadyState:
      case oExpectation:
        return cost;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms,
                                   bool is_matlab) const
{
  NodeID this2 = const_cast<UnaryOpNode *>(this);

  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      arg->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
        temporary_terms.insert(this2);
    }
}

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms,
                                   map<NodeID, pair<int, int> > &first_occurence,
                                   int Curr_block,
                                   vector< vector<temporary_terms_type> > &v_temporary_terms,
                                   int equation) const
{
  NodeID this2 = const_cast<UnaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

void
UnaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    arg->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
}

void
UnaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                         const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  // Always put parenthesis around uminus nodes
  if (op_code == oUminus)
    output << LEFT_PAR(output_type);

  switch (op_code)
    {
    case oUminus:
      output << "-";
      break;
    case oExp:
      output << "exp";
      break;
    case oLog:
      output << "log";
      break;
    case oLog10:
      if (IS_LATEX(output_type))
        output << "log_{10}";
      else
        output << "log10";
      break;
    case oCos:
      output << "cos";
      break;
    case oSin:
      output << "sin";
      break;
    case oTan:
      output << "tan";
      break;
    case oAcos:
      output << "acos";
      break;
    case oAsin:
      output << "asin";
      break;
    case oAtan:
      output << "atan";
      break;
    case oCosh:
      output << "cosh";
      break;
    case oSinh:
      output << "sinh";
      break;
    case oTanh:
      output << "tanh";
      break;
    case oAcosh:
      output << "acosh";
      break;
    case oAsinh:
      output << "asinh";
      break;
    case oAtanh:
      output << "atanh";
      break;
    case oSqrt:
      output << "sqrt";
      break;
    case oSteadyState:
      ExprNodeOutputType new_output_type;
      switch (output_type)
        {
        case oMatlabDynamicModel:
          new_output_type = oMatlabDynamicSteadyStateOperator;
          break;
        case oLatexDynamicModel:
          new_output_type = oLatexDynamicSteadyStateOperator;
          break;
        case oCDynamicModel:
          new_output_type = oCDynamicSteadyStateOperator;
          break;
        case oMatlabDynamicModelSparse:
        case oMatlabDynamicModelSparseLocalTemporaryTerms:
          new_output_type = oMatlabDynamicSparseSteadyStateOperator;
          break;
        default:
          new_output_type = output_type;
          break;
        }
      arg->writeOutput(output, new_output_type, temporary_terms);
      return;
    case oExpectation:
      assert(0);
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != oUminus
      || (op_code == oUminus
          && arg->precedence(output_type, temporary_terms) < precedence(output_type, temporary_terms)))
    {
      output << LEFT_PAR(output_type);
      close_parenthesis = true;
    }

  // Write argument
  arg->writeOutput(output, output_type, temporary_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  // Close parenthesis for uminus
  if (op_code == oUminus)
    output << RIGHT_PAR(output_type);
}

double
UnaryOpNode::eval_opcode(UnaryOpcode op_code, double v) throw (EvalException)
{
  switch (op_code)
    {
    case oUminus:
      return (-v);
    case oExp:
      return (exp(v));
    case oLog:
      return (log(v));
    case oLog10:
      return (log10(v));
    case oCos:
      return (cos(v));
    case oSin:
      return (sin(v));
    case oTan:
      return (tan(v));
    case oAcos:
      return (acos(v));
    case oAsin:
      return (asin(v));
    case oAtan:
      return (atan(v));
    case oCosh:
      return (cosh(v));
    case oSinh:
      return (sinh(v));
    case oTanh:
      return (tanh(v));
    case oAcosh:
      return (acosh(v));
    case oAsinh:
      return (asinh(v));
    case oAtanh:
      return (atanh(v));
    case oSqrt:
      return (sqrt(v));
    case oSteadyState:
      return (v);
    case oExpectation:
      throw EvalException();
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
UnaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v = arg->eval(eval_context);

  return eval_opcode(op_code, v);
}

void
UnaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(map_idx[idx]);
          fldt.write(CompileCode);
        }
      else
        {
          FLDST_ fldst(map_idx[idx]);
          fldst.write(CompileCode);
        }
      return;
    }
  if (op_code == oSteadyState)
    arg->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, true);
  else
    {
      arg->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
      FUNARY_ funary(op_code);
      funary.write(CompileCode);
    }
}

void
UnaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg->collectVariables(type_arg, result);
}

pair<int, NodeID>
UnaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  pair<bool, NodeID > res = arg->normalizeEquation(var_endo, List_of_Op_RHS);
  int is_endogenous_present = res.first;
  NodeID New_NodeID = res.second;
  /*if(res.second.second)*/
  if (is_endogenous_present == 2)
    return (make_pair(2, (NodeID) NULL));
  else if (is_endogenous_present)
    {
      switch (op_code)
        {
        case oUminus:
          List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((NodeID) NULL, (NodeID) NULL)));
          return (make_pair(1, (NodeID) NULL));
        case oExp:
          List_of_Op_RHS.push_back(make_pair(oLog, make_pair((NodeID) NULL, (NodeID) NULL)));
          return (make_pair(1, (NodeID) NULL));
        case oLog:
          List_of_Op_RHS.push_back(make_pair(oExp, make_pair((NodeID) NULL, (NodeID) NULL)));
          return (make_pair(1, (NodeID) NULL));
        case oLog10:
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair((NodeID) NULL, datatree.AddNumConstant("10"))));
          return (make_pair(1, (NodeID) NULL));
        case oCos:
          return (make_pair(1, (NodeID) NULL));
        case oSin:
          return (make_pair(1, (NodeID) NULL));
        case oTan:
          return (make_pair(1, (NodeID) NULL));
        case oAcos:
          return (make_pair(1, (NodeID) NULL));
        case oAsin:
          return (make_pair(1, (NodeID) NULL));
        case oAtan:
          return (make_pair(1, (NodeID) NULL));
        case oCosh:
          return (make_pair(1, (NodeID) NULL));
        case oSinh:
          return (make_pair(1, (NodeID) NULL));
        case oTanh:
          return (make_pair(1, (NodeID) NULL));
        case oAcosh:
          return (make_pair(1, (NodeID) NULL));
        case oAsinh:
          return (make_pair(1, (NodeID) NULL));
        case oAtanh:
          return (make_pair(1, (NodeID) NULL));
        case oSqrt:
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair((NodeID) NULL, datatree.AddNumConstant("2"))));
          return (make_pair(1, (NodeID) NULL));
        case oSteadyState:
          return (make_pair(1, (NodeID) NULL));
        case oExpectation:
          assert(0);
        }
    }
  else
    {
      switch (op_code)
        {
        case oUminus:
          return (make_pair(0, datatree.AddUMinus(New_NodeID)));
        case oExp:
          return (make_pair(0, datatree.AddExp(New_NodeID)));
        case oLog:
          return (make_pair(0, datatree.AddLog(New_NodeID)));
        case oLog10:
          return (make_pair(0, datatree.AddLog10(New_NodeID)));
        case oCos:
          return (make_pair(0, datatree.AddCos(New_NodeID)));
        case oSin:
          return (make_pair(0, datatree.AddSin(New_NodeID)));
        case oTan:
          return (make_pair(0, datatree.AddTan(New_NodeID)));
        case oAcos:
          return (make_pair(0, datatree.AddAcos(New_NodeID)));
        case oAsin:
          return (make_pair(0, datatree.AddAsin(New_NodeID)));
        case oAtan:
          return (make_pair(0, datatree.AddAtan(New_NodeID)));
        case oCosh:
          return (make_pair(0, datatree.AddCosh(New_NodeID)));
        case oSinh:
          return (make_pair(0, datatree.AddSinh(New_NodeID)));
        case oTanh:
          return (make_pair(0, datatree.AddTanh(New_NodeID)));
        case oAcosh:
          return (make_pair(0, datatree.AddAcosh(New_NodeID)));
        case oAsinh:
          return (make_pair(0, datatree.AddAsinh(New_NodeID)));
        case oAtanh:
          return (make_pair(0, datatree.AddAtanh(New_NodeID)));
        case oSqrt:
          return (make_pair(0, datatree.AddSqrt(New_NodeID)));
        case oSteadyState:
          return (make_pair(0, datatree.AddSteadyState(New_NodeID)));
        case oExpectation:
          assert(0);
        }
    }
  return (make_pair(1, (NodeID) NULL));
}

NodeID
UnaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg = arg->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg);
}

NodeID
UnaryOpNode::buildSimilarUnaryOpNode(NodeID alt_arg, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case oUminus:
      return alt_datatree.AddUMinus(alt_arg);
    case oExp:
      return alt_datatree.AddExp(alt_arg);
    case oLog:
      return alt_datatree.AddLog(alt_arg);
    case oLog10:
      return alt_datatree.AddLog10(alt_arg);
    case oCos:
      return alt_datatree.AddCos(alt_arg);
    case oSin:
      return alt_datatree.AddSin(alt_arg);
    case oTan:
      return alt_datatree.AddTan(alt_arg);
    case oAcos:
      return alt_datatree.AddAcos(alt_arg);
    case oAsin:
      return alt_datatree.AddAsin(alt_arg);
    case oAtan:
      return alt_datatree.AddAtan(alt_arg);
    case oCosh:
      return alt_datatree.AddCosh(alt_arg);
    case oSinh:
      return alt_datatree.AddSinh(alt_arg);
    case oTanh:
      return alt_datatree.AddTanh(alt_arg);
    case oAcosh:
      return alt_datatree.AddAcosh(alt_arg);
    case oAsinh:
      return alt_datatree.AddAsinh(alt_arg);
    case oAtanh:
      return alt_datatree.AddAtanh(alt_arg);
    case oSqrt:
      return alt_datatree.AddSqrt(alt_arg);
    case oSteadyState:
      return alt_datatree.AddSteadyState(alt_arg);
    case oExpectation:
      return alt_datatree.AddExpectation(expectation_information_set, alt_arg);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
UnaryOpNode::toStatic(DataTree &static_datatree) const
{
  NodeID sarg = arg->toStatic(static_datatree);
  return buildSimilarUnaryOpNode(sarg, static_datatree);
}

int
UnaryOpNode::maxEndoLead() const
{
  return arg->maxEndoLead();
}

int
UnaryOpNode::maxExoLead() const
{
  return arg->maxExoLead();
}

int
UnaryOpNode::maxEndoLag() const
{
  return arg->maxEndoLag();
}

int
UnaryOpNode::maxExoLag() const
{
  return arg->maxExoLag();
}

NodeID
UnaryOpNode::decreaseLeadsLags(int n) const
{
  NodeID argsubst = arg->decreaseLeadsLags(n);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

NodeID
UnaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  NodeID argsubst = arg->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

NodeID
UnaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (op_code == oUminus)
    {
      NodeID argsubst = arg->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
  else
    {
      if (maxEndoLead() >= 2)
        return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      else
        return const_cast<UnaryOpNode *>(this);
    }
}

NodeID
UnaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID argsubst = arg->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

NodeID
UnaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (op_code == oUminus)
    {
      NodeID argsubst = arg->substituteExoLead(subst_table, neweqs);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
  else
    {
      if (maxExoLead() >= 1)
        return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
      else
        return const_cast<UnaryOpNode *>(this);
    }
}

NodeID
UnaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID argsubst = arg->substituteExoLag(subst_table, neweqs);
  return buildSimilarUnaryOpNode(argsubst, datatree);
}

NodeID
UnaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  switch (op_code)
    {
    case oExpectation:
      {
        subst_table_t::iterator it = subst_table.find(const_cast<UnaryOpNode *>(this));

        if (it != subst_table.end())
          return const_cast<VariableNode *>(it->second);

        //Arriving here, we need to create an auxiliary variable for this Expectation Operator:
        int symb_id = datatree.symbol_table.addExpectationAuxiliaryVar(expectation_information_set, arg->idx); //AUXE_period_arg.idx
        NodeID newAuxE = datatree.AddVariable(symb_id, 0);

        if (partial_information_model && expectation_information_set == 0)
          {
            //Ensure x is a single variable as opposed to an expression
            if (dynamic_cast<VariableNode *>(arg) == NULL)
              {
                cerr << "ERROR: In Partial Information models, EXPECTATION(0)(X) can only be used when X is a single variable." << endl;
                exit(EXIT_FAILURE);
              }
          }
        else
          {
            //take care of any nested expectation operators by calling arg->substituteExpectation(.), then decreaseLeadsLags for this oExpectation operator
            //arg(lag-period) (holds entire subtree of arg(lag-period)
            NodeID substexpr = (arg->substituteExpectation(subst_table, neweqs, partial_information_model))->decreaseLeadsLags(expectation_information_set);
            assert(substexpr != NULL);

            neweqs.push_back(dynamic_cast<BinaryOpNode *>(datatree.AddEqual(newAuxE, substexpr))); //AUXE_period_arg.idx = arg(lag-period)

            newAuxE = datatree.AddVariable(symb_id, expectation_information_set);
          }

        assert(dynamic_cast<VariableNode *>(newAuxE) != NULL);
        subst_table[this] = dynamic_cast<VariableNode *>(newAuxE);
        return newAuxE;
      }
    default:
      NodeID argsubst = arg->substituteExpectation(subst_table, neweqs, partial_information_model);
      return buildSimilarUnaryOpNode(argsubst, datatree);
    }
}

bool
UnaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
UnaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
                           BinaryOpcode op_code_arg, const NodeID arg2_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  op_code(op_code_arg)
{
  datatree.binary_op_node_map[make_pair(make_pair(arg1, arg2), op_code)] = this;
}

void
BinaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg1->prepareForDerivation();
  arg2->prepareForDerivation();

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg1->non_null_derivatives and arg2->non_null_derivatives
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

NodeID
BinaryOpNode::composeDerivatives(NodeID darg1, NodeID darg2)
{
  NodeID t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case oPlus:
      return datatree.AddPlus(darg1, darg2);
    case oMinus:
      return datatree.AddMinus(darg1, darg2);
    case oTimes:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      return datatree.AddPlus(t11, t12);
    case oDivide:
      if (darg2 != datatree.Zero)
        {
          t11 = datatree.AddTimes(darg1, arg2);
          t12 = datatree.AddTimes(darg2, arg1);
          t13 = datatree.AddMinus(t11, t12);
          t14 = datatree.AddTimes(arg2, arg2);
          return datatree.AddDivide(t13, t14);
        }
      else
        return datatree.AddDivide(darg1, arg2);
    case oLess:
    case oGreater:
    case oLessEqual:
    case oGreaterEqual:
    case oEqualEqual:
    case oDifferent:
      return datatree.Zero;
    case oPower:
      if (darg2 == datatree.Zero)
        {
          if (darg1 == datatree.Zero)
            return datatree.Zero;
          else
            {
              t11 = datatree.AddMinus(arg2, datatree.One);
              t12 = datatree.AddPower(arg1, t11);
              t13 = datatree.AddTimes(arg2, t12);
              return datatree.AddTimes(darg1, t13);
            }
        }
      else
        {
          t11 = datatree.AddLog(arg1);
          t12 = datatree.AddTimes(darg2, t11);
          t13 = datatree.AddTimes(darg1, arg2);
          t14 = datatree.AddDivide(t13, arg1);
          t15 = datatree.AddPlus(t12, t14);
          return datatree.AddTimes(t15, this);
        }
    case oMax:
      t11 = datatree.AddGreater(arg1, arg2);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    case oMin:
      t11 = datatree.AddGreater(arg2, arg1);
      t12 = datatree.AddTimes(t11, darg1);
      t13 = datatree.AddMinus(datatree.One, t11);
      t14 = datatree.AddTimes(t13, darg2);
      return datatree.AddPlus(t14, t12);
    case oEqual:
      return datatree.AddMinus(darg1, darg2);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
BinaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg1 = arg1->getDerivative(deriv_id);
  NodeID darg2 = arg2->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2);
}

int
BinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case oEqual:
      return 0;
    case oEqualEqual:
    case oDifferent:
      return 1;
    case oLessEqual:
    case oGreaterEqual:
    case oLess:
    case oGreater:
      return 2;
    case oPlus:
    case oMinus:
      return 3;
    case oTimes:
    case oDivide:
      return 4;
    case oPower:
      if (IS_C(output_type))
        // In C, power operator is of the form pow(a, b)
        return 100;
      else
        return 5;
    case oMin:
    case oMax:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
BinaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // For a temporary term, the cost is null
  if (it != temporary_terms.end())
    return 0;

  int cost = arg1->cost(temporary_terms, is_matlab);
  cost += arg2->cost(temporary_terms, is_matlab);

  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oLess:
      case oGreater:
      case oLessEqual:
      case oGreaterEqual:
      case oEqualEqual:
      case oDifferent:
        return cost + 60;
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 90;
      case oMax:
      case oMin:
        return cost + 110;
      case oDivide:
        return cost + 990;
      case oPower:
        return cost + 1160;
      case oEqual:
        return cost;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oLess:
      case oGreater:
      case oLessEqual:
      case oGreaterEqual:
      case oEqualEqual:
      case oDifferent:
        return cost + 2;
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 4;
      case oMax:
      case oMin:
        return cost + 5;
      case oDivide:
        return cost + 15;
      case oPower:
        return cost + 520;
      case oEqual:
        return cost;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    bool is_matlab) const
{
  NodeID this2 = const_cast<BinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = 1;
      arg1->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
    }
  else
    {
      /* If the node has already been encountered, increment its ref count
         and declare it as a temporary term if it is too costly (except if it is
         an equal node: we don't want them as temporary terms) */
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab)
          && op_code != oEqual)
        temporary_terms.insert(this2);
    }
}

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    map<NodeID, pair<int, int> > &first_occurence,
                                    int Curr_block,
                                    vector<vector<temporary_terms_type> > &v_temporary_terms,
                                    int equation) const
{
  NodeID this2 = const_cast<BinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C
          && op_code != oEqual)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

double
BinaryOpNode::eval_opcode(double v1, BinaryOpcode op_code, double v2) throw (EvalException)
{
  switch (op_code)
    {
    case oPlus:
      return (v1 + v2);
    case oMinus:
      return (v1 - v2);
    case oTimes:
      return (v1 * v2);
    case oDivide:
      return (v1 / v2);
    case oPower:
      return (pow(v1, v2));
    case oMax:
      if (v1 < v2)
        return v2;
      else
        return v1;
    case oMin:
      if (v1 > v2)
        return v2;
      else
        return v1;
    case oLess:
      return (v1 < v2);
    case oGreater:
      return (v1 > v2);
    case oLessEqual:
      return (v1 <= v2);
    case oGreaterEqual:
      return (v1 >= v2);
    case oEqualEqual:
      return (v1 == v2);
    case oDifferent:
      return (v1 != v2);
    case oEqual:
      throw EvalException();
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
BinaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);

  return eval_opcode(v1, op_code, v2);
}

void
BinaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(map_idx[idx]);
          fldt.write(CompileCode);
        }
      else
        {
          FLDST_ fldst(map_idx[idx]);
          fldst.write(CompileCode);
        }
      return;
    }
  arg1->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  arg2->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  FBINARY_ fbinary(op_code);
  fbinary.write(CompileCode);
}

void
BinaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      arg1->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg2->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
    }
}

void
BinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type == oMatlabDynamicModelSparse)
        output << "T" << idx << "(it_)";
      else
        output << "T" << idx;
      return;
    }

  // Treat special case of power operator in C, and case of max and min operators
  if ((op_code == oPower && IS_C(output_type)) || op_code == oMax || op_code == oMin)
    {
      switch (op_code)
        {
        case oPower:
          output << "pow(";
          break;
        case oMax:
          output << "max(";
          break;
        case oMin:
          output << "min(";
          break;
        default:
          ;
        }
      arg1->writeOutput(output, output_type, temporary_terms);
      output << ",";
      arg2->writeOutput(output, output_type, temporary_terms);
      output << ")";
      return;
    }

  int prec = precedence(output_type, temporary_terms);

  bool close_parenthesis = false;

  if (IS_LATEX(output_type) && op_code == oDivide)
    output << "\\frac{";
  else
    {
      // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
      BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
      if (arg1->precedence(output_type, temporary_terms) < prec
          || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write left argument
  arg1->writeOutput(output, output_type, temporary_terms);

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);

  if (IS_LATEX(output_type) && op_code == oDivide)
    output << "}";

  // Write current operator symbol
  switch (op_code)
    {
    case oPlus:
      output << "+";
      break;
    case oMinus:
      output << "-";
      break;
    case oTimes:
      if (IS_LATEX(output_type))
        output << "\\, ";
      else
        output << "*";
      break;
    case oDivide:
      if (!IS_LATEX(output_type))
        output << "/";
      break;
    case oPower:
      output << "^";
      break;
    case oLess:
      output << "<";
      break;
    case oGreater:
      output << ">";
      break;
    case oLessEqual:
      if (IS_LATEX(output_type))
        output << "\\leq ";
      else
        output << "<=";
      break;
    case oGreaterEqual:
      if (IS_LATEX(output_type))
        output << "\\geq ";
      else
        output << ">=";
      break;
    case oEqualEqual:
      output << "==";
      break;
    case oDifferent:
      if (IS_MATLAB(output_type))
        output << "~=";
      else
        {
          if (IS_C(output_type))
            output << "!=";
          else
            output << "\\neq ";
        }
      break;
    case oEqual:
      output << "=";
      break;
    default:
      ;
    }

  close_parenthesis = false;

  if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
    output << "{";
  else
    {
      /* Add parenthesis around right argument if:
         - its precedence is lower than those of the current node
         - it is a power operator and current operator is also a power operator
         - it is a minus operator with same precedence than current operator
         - it is a divide operator with same precedence than current operator */
      BinaryOpNode *barg2 = dynamic_cast<BinaryOpNode *>(arg2);
      int arg2_prec = arg2->precedence(output_type, temporary_terms);
      if (arg2_prec < prec
          || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower && !IS_LATEX(output_type))
          || (op_code == oMinus && arg2_prec == prec)
          || (op_code == oDivide && arg2_prec == prec && !IS_LATEX(output_type)))
        {
          output << LEFT_PAR(output_type);
          close_parenthesis = true;
        }
    }

  // Write right argument
  arg2->writeOutput(output, output_type, temporary_terms);

  if (IS_LATEX(output_type) && (op_code == oPower || op_code == oDivide))
    output << "}";

  if (close_parenthesis)
    output << RIGHT_PAR(output_type);
}

void
BinaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectVariables(type_arg, result);
  arg2->collectVariables(type_arg, result);
}

NodeID
BinaryOpNode::Compute_RHS(NodeID arg1, NodeID arg2, int op, int op_type) const
{
  temporary_terms_type temp;
  switch (op_type)
    {
    case 0: /*Unary Operator*/
      switch (op)
        {
        case oUminus:
          return (datatree.AddUMinus(arg1));
          break;
        case oExp:
          return (datatree.AddExp(arg1));
          break;
        case oLog:
          return (datatree.AddLog(arg1));
          break;
        case oLog10:
          return (datatree.AddLog10(arg1));
          break;
        }
      break;
    case 1: /*Binary Operator*/
      switch (op)
        {
        case oPlus:
          return (datatree.AddPlus(arg1, arg2));
          break;
        case oMinus:
          return (datatree.AddMinus(arg1, arg2));
          break;
        case oTimes:
          return (datatree.AddTimes(arg1, arg2));
          break;
        case oDivide:
          return (datatree.AddDivide(arg1, arg2));
          break;
        case oPower:
          return (datatree.AddPower(arg1, arg2));
          break;
        }
      break;
    }
  return ((NodeID) NULL);
}

pair<int, NodeID>
BinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  vector<pair<int, pair<NodeID, NodeID> > > List_of_Op_RHS1, List_of_Op_RHS2;
  int is_endogenous_present_1, is_endogenous_present_2;
  pair<int, NodeID> res;
  NodeID NodeID_1, NodeID_2;
  res = arg1->normalizeEquation(var_endo, List_of_Op_RHS1);
  is_endogenous_present_1 = res.first;
  NodeID_1 = res.second;

  res = arg2->normalizeEquation(var_endo, List_of_Op_RHS2);
  is_endogenous_present_2 = res.first;
  NodeID_2 = res.second;
  if (is_endogenous_present_1 == 2 || is_endogenous_present_2 == 2)
    return (make_pair(2, (NodeID) NULL));
  else if (is_endogenous_present_1 && is_endogenous_present_2)
    return (make_pair(2, (NodeID) NULL));
  else if (is_endogenous_present_1)
    {
      if (op_code == oEqual)
        {
          pair<int, pair<NodeID, NodeID> > it;
          int oo = List_of_Op_RHS1.size();
          for (int i = 0; i < oo; i++)
            {
              it = List_of_Op_RHS1.back();
              List_of_Op_RHS1.pop_back();
              if (it.second.first && !it.second.second) /*Binary operator*/
                NodeID_2 = Compute_RHS(NodeID_2, (BinaryOpNode *) it.second.first, it.first, 1);
              else if (it.second.second && !it.second.first) /*Binary operator*/
                NodeID_2 = Compute_RHS(it.second.second, NodeID_2, it.first, 1);
              else if (it.second.second && it.second.first) /*Binary operator*/
                NodeID_2 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
              else                                          /*Unary operator*/
                NodeID_2 = Compute_RHS((UnaryOpNode *) NodeID_2, (UnaryOpNode *) it.second.first, it.first, 0);
            }
        }
      else
        List_of_Op_RHS = List_of_Op_RHS1;
    }
  else if (is_endogenous_present_2)
    {
      if (op_code == oEqual)
        {
          int oo = List_of_Op_RHS2.size();
          for (int i = 0; i < oo; i++)
            {
              pair<int, pair<NodeID, NodeID> > it;
              it = List_of_Op_RHS2.back();
              List_of_Op_RHS2.pop_back();
              if (it.second.first && !it.second.second) /*Binary operator*/
                NodeID_1 = Compute_RHS((BinaryOpNode *) NodeID_1, (BinaryOpNode *) it.second.first, it.first, 1);
              else if (it.second.second && !it.second.first) /*Binary operator*/
                NodeID_1 = Compute_RHS((BinaryOpNode *) it.second.second, (BinaryOpNode *) NodeID_1, it.first, 1);
              else if (it.second.second && it.second.first) /*Binary operator*/
                NodeID_1 = Compute_RHS(it.second.first, it.second.second, it.first, 1);
              else
                NodeID_1 = Compute_RHS((UnaryOpNode *) NodeID_1, (UnaryOpNode *) it.second.first, it.first, 0);
            }
        }
      else
        List_of_Op_RHS = List_of_Op_RHS2;
    }
  switch (op_code)
    {
    case oPlus:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddPlus(NodeID_1, NodeID_2), (NodeID) NULL)));
          return (make_pair(0, datatree.AddPlus(NodeID_1, NodeID_2)));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        return (make_pair(1, (NodeID) NULL));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_1, (NodeID) NULL)));
          return (make_pair(1, NodeID_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_2, (NodeID) NULL)));
          return (make_pair(1, NodeID_2));
        }
      break;
    case oMinus:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(datatree.AddMinus(NodeID_1, NodeID_2), (NodeID) NULL)));
          return (make_pair(0, datatree.AddMinus(NodeID_1, NodeID_2)));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        return (make_pair(1, (NodeID) NULL));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oUminus, make_pair((NodeID) NULL, (NodeID) NULL)));
          List_of_Op_RHS.push_back(make_pair(oMinus, make_pair(NodeID_1, (NodeID) NULL)));
          return (make_pair(1, NodeID_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oPlus, make_pair(NodeID_2, (NodeID) NULL)));
          return (make_pair(1, datatree.AddUMinus(NodeID_2)));
        }
      break;
    case oTimes:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddTimes(NodeID_1, NodeID_2)));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(NodeID_1, (NodeID) NULL)));
          return (make_pair(1, NodeID_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair(NodeID_2, (NodeID) NULL)));
          return (make_pair(1, NodeID_2));
        }
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oDivide:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddDivide(NodeID_1, NodeID_2)));
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oDivide, make_pair((NodeID) NULL, NodeID_1)));
          return (make_pair(1, NodeID_1));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oTimes, make_pair(NodeID_2, (NodeID) NULL)));
          return (make_pair(1, NodeID_2));
        }
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oPower:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddPower(NodeID_1, NodeID_2)));
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          List_of_Op_RHS.push_back(make_pair(oPower, make_pair(datatree.AddDivide(datatree.AddNumConstant("1"), NodeID_2), (NodeID) NULL)));
          return (make_pair(1, (NodeID) NULL));
        }
      break;
    case oEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), datatree.AddMinus(NodeID_2, NodeID_1))
                            ));
        }
      else if (is_endogenous_present_1 && is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), datatree.Zero)
                            ));
        }
      else if (!is_endogenous_present_1 && is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), /*datatree.AddUMinus(NodeID_1)*/ NodeID_1)
                            ));
        }
      else if (is_endogenous_present_1 && !is_endogenous_present_2)
        {
          return (make_pair(0,
                            datatree.AddEqual(datatree.AddVariable(datatree.symbol_table.getName(datatree.symbol_table.getID(eEndogenous, var_endo)), 0), NodeID_2)
                            ));
        }
      break;
    case oMax:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddMax(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oMin:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddMin(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oLess:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddLess(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oGreater:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddGreater(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oLessEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddLessEqual(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oGreaterEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddGreaterEqual(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oEqualEqual:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddEqualEqual(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    case oDifferent:
      if (!is_endogenous_present_1 && !is_endogenous_present_2)
        return (make_pair(0, datatree.AddDifferent(NodeID_1, NodeID_2)));
      else
        return (make_pair(1, (NodeID) NULL));
      break;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
BinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2);
}

NodeID
BinaryOpNode::buildSimilarBinaryOpNode(NodeID alt_arg1, NodeID alt_arg2, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case oPlus:
      return alt_datatree.AddPlus(alt_arg1, alt_arg2);
    case oMinus:
      return alt_datatree.AddMinus(alt_arg1, alt_arg2);
    case oTimes:
      return alt_datatree.AddTimes(alt_arg1, alt_arg2);
    case oDivide:
      return alt_datatree.AddDivide(alt_arg1, alt_arg2);
    case oPower:
      return alt_datatree.AddPower(alt_arg1, alt_arg2);
    case oEqual:
      return alt_datatree.AddEqual(alt_arg1, alt_arg2);
    case oMax:
      return alt_datatree.AddMax(alt_arg1, alt_arg2);
    case oMin:
      return alt_datatree.AddMin(alt_arg1, alt_arg2);
    case oLess:
      return alt_datatree.AddLess(alt_arg1, alt_arg2);
    case oGreater:
      return alt_datatree.AddGreater(alt_arg1, alt_arg2);
    case oLessEqual:
      return alt_datatree.AddLessEqual(alt_arg1, alt_arg2);
    case oGreaterEqual:
      return alt_datatree.AddGreaterEqual(alt_arg1, alt_arg2);
    case oEqualEqual:
      return alt_datatree.AddEqualEqual(alt_arg1, alt_arg2);
    case oDifferent:
      return alt_datatree.AddDifferent(alt_arg1, alt_arg2);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
BinaryOpNode::toStatic(DataTree &static_datatree) const
{
  NodeID sarg1 = arg1->toStatic(static_datatree);
  NodeID sarg2 = arg2->toStatic(static_datatree);
  return buildSimilarBinaryOpNode(sarg1, sarg2, static_datatree);
}

int
BinaryOpNode::maxEndoLead() const
{
  return max(arg1->maxEndoLead(), arg2->maxEndoLead());
}

int
BinaryOpNode::maxExoLead() const
{
  return max(arg1->maxExoLead(), arg2->maxExoLead());
}

int
BinaryOpNode::maxEndoLag() const
{
  return max(arg1->maxEndoLag(), arg2->maxEndoLag());
}

int
BinaryOpNode::maxExoLag() const
{
  return max(arg1->maxExoLag(), arg2->maxExoLag());
}

NodeID
BinaryOpNode::decreaseLeadsLags(int n) const
{
  NodeID arg1subst = arg1->decreaseLeadsLags(n);
  NodeID arg2subst = arg2->decreaseLeadsLags(n);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

NodeID
BinaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  NodeID arg1subst = arg1->decreaseLeadsLagsPredeterminedVariables();
  NodeID arg2subst = arg2->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

NodeID
BinaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst, arg2subst;
  int maxendolead1 = arg1->maxEndoLead(), maxendolead2 = arg2->maxEndoLead();

  if (maxendolead1 < 2 && maxendolead2 < 2)
    return const_cast<BinaryOpNode *>(this);

  switch (op_code)
    {
    case oPlus:
    case oMinus:
    case oEqual:
      arg1subst = maxendolead1 >= 2 ? arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs) : arg1;
      arg2subst = maxendolead2 >= 2 ? arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs) : arg2;
      return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
    case oTimes:
    case oDivide:
      if (maxendolead1 >= 2 && maxendolead2 == 0 && arg2->maxExoLead() == 0)
        {
          arg1subst = arg1->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
        }
      if (maxendolead1 == 0 && arg1->maxExoLead() == 0
          && maxendolead2 >= 2 && op_code == oTimes)
        {
          arg2subst = arg2->substituteEndoLeadGreaterThanTwo(subst_table, neweqs);
          return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
        }
      return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    default:
      return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    }
}

NodeID
BinaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst = arg1->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  NodeID arg2subst = arg2->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

NodeID
BinaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst, arg2subst;
  int maxexolead1 = arg1->maxExoLead(), maxexolead2 = arg2->maxExoLead();

  if (maxexolead1 < 1 && maxexolead2 < 1)
    return const_cast<BinaryOpNode *>(this);

  switch (op_code)
    {
    case oPlus:
    case oMinus:
    case oEqual:
      arg1subst = maxexolead1 >= 1 ? arg1->substituteExoLead(subst_table, neweqs) : arg1;
      arg2subst = maxexolead2 >= 1 ? arg2->substituteExoLead(subst_table, neweqs) : arg2;
      return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
    case oTimes:
    case oDivide:
      if (maxexolead1 >= 1 && maxexolead2 == 0 && arg2->maxEndoLead() == 0)
        {
          arg1subst = arg1->substituteExoLead(subst_table, neweqs);
          return buildSimilarBinaryOpNode(arg1subst, arg2, datatree);
        }
      if (maxexolead1 == 0 && arg1->maxEndoLead() == 0
          && maxexolead2 >= 1 && op_code == oTimes)
        {
          arg2subst = arg2->substituteExoLead(subst_table, neweqs);
          return buildSimilarBinaryOpNode(arg1, arg2subst, datatree);
        }
      return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    default:
      return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
    }
}

NodeID
BinaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst = arg1->substituteExoLag(subst_table, neweqs);
  NodeID arg2subst = arg2->substituteExoLag(subst_table, neweqs);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

NodeID
BinaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  NodeID arg1subst = arg1->substituteExpectation(subst_table, neweqs, partial_information_model);
  NodeID arg2subst = arg2->substituteExpectation(subst_table, neweqs, partial_information_model);
  return buildSimilarBinaryOpNode(arg1subst, arg2subst, datatree);
}

bool
BinaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
BinaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

TrinaryOpNode::TrinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
                             TrinaryOpcode op_code_arg, const NodeID arg2_arg, const NodeID arg3_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  arg3(arg3_arg),
  op_code(op_code_arg)
{
  datatree.trinary_op_node_map[make_pair(make_pair(make_pair(arg1, arg2), arg3), op_code)] = this;
}

void
TrinaryOpNode::prepareForDerivation()
{
  if (preparedForDerivation)
    return;

  preparedForDerivation = true;

  arg1->prepareForDerivation();
  arg2->prepareForDerivation();
  arg3->prepareForDerivation();

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg{1,2,3}->non_null_derivatives
  set<int> non_null_derivatives_tmp;
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives_tmp, non_null_derivatives_tmp.begin()));
  set_union(non_null_derivatives_tmp.begin(),
            non_null_derivatives_tmp.end(),
            arg3->non_null_derivatives.begin(),
            arg3->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

NodeID
TrinaryOpNode::composeDerivatives(NodeID darg1, NodeID darg2, NodeID darg3)
{

  NodeID t11, t12, t13, t14, t15;

  switch (op_code)
    {
    case oNormcdf:
      // normal pdf is inlined in the tree
      NodeID y;
      // sqrt(2*pi)
      t14 = datatree.AddSqrt(datatree.AddTimes(datatree.Two, datatree.Pi));
      // x - mu
      t12 = datatree.AddMinus(arg1, arg2);
      // y = (x-mu)/sigma
      y = datatree.AddDivide(t12, arg3);
      // (x-mu)^2/sigma^2
      t12 = datatree.AddTimes(y, y);
      // -(x-mu)^2/sigma^2
      t13 = datatree.AddUMinus(t12);
      // -((x-mu)^2/sigma^2)/2
      t12 = datatree.AddDivide(t13, datatree.Two);
      // exp(-((x-mu)^2/sigma^2)/2)
      t13 = datatree.AddExp(t12);
      // derivative of a standardized normal
      // t15 = (1/sqrt(2*pi))*exp(-y^2/2)
      t15 = datatree.AddDivide(t13, t14);
      // derivatives thru x
      t11 = datatree.AddDivide(darg1, arg3);
      // derivatives thru mu
      t12 = datatree.AddDivide(darg2, arg3);
      // intermediary sum
      t14 = datatree.AddMinus(t11, t12);
      // derivatives thru sigma
      t11 = datatree.AddDivide(y, arg3);
      t12 = datatree.AddTimes(t11, darg3);
      //intermediary sum
      t11 = datatree.AddMinus(t14, t12);
      // total derivative:
      // (darg1/sigma - darg2/sigma - darg3*(x-mu)/sigma^2) * t15
      // where t15 is the derivative of a standardized normal
      return datatree.AddTimes(t11, t15);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
TrinaryOpNode::computeDerivative(int deriv_id)
{
  NodeID darg1 = arg1->getDerivative(deriv_id);
  NodeID darg2 = arg2->getDerivative(deriv_id);
  NodeID darg3 = arg3->getDerivative(deriv_id);
  return composeDerivatives(darg1, darg2, darg3);
}

int
TrinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch (op_code)
    {
    case oNormcdf:
      return 100;
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

int
TrinaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  // For a temporary term, the cost is null
  if (it != temporary_terms.end())
    return 0;

  int cost = arg1->cost(temporary_terms, is_matlab);
  cost += arg2->cost(temporary_terms, is_matlab);

  if (is_matlab)
    // Cost for Matlab files
    switch (op_code)
      {
      case oNormcdf:
        return cost+1000;
      }
  else
    // Cost for C files
    switch (op_code)
      {
      case oNormcdf:
        return cost+1000;
      }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

void
TrinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     bool is_matlab) const
{
  NodeID this2 = const_cast<TrinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = 1;
      arg1->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      arg3->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
    }
  else
    {
      // If the node has already been encountered, increment its ref count
      //  and declare it as a temporary term if it is too costly
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
        temporary_terms.insert(this2);
    }
}

void
TrinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector<vector<temporary_terms_type> > &v_temporary_terms,
                                     int equation) const
{
  NodeID this2 = const_cast<TrinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = make_pair(Curr_block, equation);
      arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
      arg3->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, v_temporary_terms, equation);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          v_temporary_terms[first_occurence[this2].first][first_occurence[this2].second].insert(this2);
        }
    }
}

double
TrinaryOpNode::eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException)
{
  switch (op_code)
    {
    case oNormcdf:
      return (0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

double
TrinaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);
  double v3 = arg3->eval(eval_context);

  return eval_opcode(v1, op_code, v2, v3);
}

void
TrinaryOpNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (dynamic)
        {
          FLDT_ fldt(map_idx[idx]);
          fldt.write(CompileCode);
        }
      else
        {
          FLDST_ fldst(map_idx[idx]);
          fldst.write(CompileCode);
        }
      return;
    }
  arg1->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  arg2->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  arg3->compile(CompileCode, lhs_rhs, temporary_terms, map_idx, dynamic, steady_dynamic);
  FTRINARY_ ftrinary(op_code);
  ftrinary.write(CompileCode);
}

void
TrinaryOpNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      arg1->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg2->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
      arg3->collectTemporary_terms(temporary_terms, temporary_terms_inuse, Curr_Block);
    }
}

void
TrinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                           const temporary_terms_type &temporary_terms) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<TrinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  switch (op_code)
    {
    case oNormcdf:
      if (IS_C(output_type))
        {
          // In C, there is no normcdf() primitive, so use erf()
          output << "(0.5*(1+erf(((";
          arg1->writeOutput(output, output_type, temporary_terms);
          output << ")-(";
          arg2->writeOutput(output, output_type, temporary_terms);
          output << "))/(";
          arg3->writeOutput(output, output_type, temporary_terms);
          output << ")/M_SQRT2)))";
        }
      else
        {
          output << "normcdf(";
          arg1->writeOutput(output, output_type, temporary_terms);
          output << ",";
          arg2->writeOutput(output, output_type, temporary_terms);
          output << ",";
          arg3->writeOutput(output, output_type, temporary_terms);
          output << ")";
        }
      break;
    }
}

void
TrinaryOpNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  arg1->collectVariables(type_arg, result);
  arg2->collectVariables(type_arg, result);
  arg3->collectVariables(type_arg, result);
}

pair<int, NodeID>
TrinaryOpNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > > &List_of_Op_RHS) const
{
  pair<int, NodeID> res = arg1->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_1 = res.first;
  NodeID NodeID_1 = res.second;
  res = arg2->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_2 = res.first;
  NodeID NodeID_2 = res.second;
  res = arg3->normalizeEquation(var_endo, List_of_Op_RHS);
  bool is_endogenous_present_3 = res.first;
  NodeID NodeID_3 = res.second;
  if (!is_endogenous_present_1 && !is_endogenous_present_2 && !is_endogenous_present_3)
    return (make_pair(0, datatree.AddNormcdf(NodeID_1, NodeID_2, NodeID_3)));
  else
    return (make_pair(1, (NodeID) NULL));
}

NodeID
TrinaryOpNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  NodeID darg1 = arg1->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg2 = arg2->getChainRuleDerivative(deriv_id, recursive_variables);
  NodeID darg3 = arg3->getChainRuleDerivative(deriv_id, recursive_variables);
  return composeDerivatives(darg1, darg2, darg3);
}

NodeID
TrinaryOpNode::buildSimilarTrinaryOpNode(NodeID alt_arg1, NodeID alt_arg2, NodeID alt_arg3, DataTree &alt_datatree) const
{
  switch (op_code)
    {
    case oNormcdf:
      return alt_datatree.AddNormcdf(alt_arg1, alt_arg2, alt_arg3);
    }
  // Suppress GCC warning
  exit(EXIT_FAILURE);
}

NodeID
TrinaryOpNode::toStatic(DataTree &static_datatree) const
{
  NodeID sarg1 = arg1->toStatic(static_datatree);
  NodeID sarg2 = arg2->toStatic(static_datatree);
  NodeID sarg3 = arg3->toStatic(static_datatree);
  return buildSimilarTrinaryOpNode(sarg1, sarg2, sarg3, static_datatree);
}

int
TrinaryOpNode::maxEndoLead() const
{
  return max(arg1->maxEndoLead(), max(arg2->maxEndoLead(), arg3->maxEndoLead()));
}

int
TrinaryOpNode::maxExoLead() const
{
  return max(arg1->maxExoLead(), max(arg2->maxExoLead(), arg3->maxExoLead()));
}

int
TrinaryOpNode::maxEndoLag() const
{
  return max(arg1->maxEndoLag(), max(arg2->maxEndoLag(), arg3->maxEndoLag()));
}

int
TrinaryOpNode::maxExoLag() const
{
  return max(arg1->maxExoLag(), max(arg2->maxExoLag(), arg3->maxExoLag()));
}

NodeID
TrinaryOpNode::decreaseLeadsLags(int n) const
{
  NodeID arg1subst = arg1->decreaseLeadsLags(n);
  NodeID arg2subst = arg2->decreaseLeadsLags(n);
  NodeID arg3subst = arg3->decreaseLeadsLags(n);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

NodeID
TrinaryOpNode::decreaseLeadsLagsPredeterminedVariables() const
{
  NodeID arg1subst = arg1->decreaseLeadsLagsPredeterminedVariables();
  NodeID arg2subst = arg2->decreaseLeadsLagsPredeterminedVariables();
  NodeID arg3subst = arg3->decreaseLeadsLagsPredeterminedVariables();
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

NodeID
TrinaryOpNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (maxEndoLead() < 2)
    return const_cast<TrinaryOpNode *>(this);
  else
    return createEndoLeadAuxiliaryVarForMyself(subst_table, neweqs);
}

NodeID
TrinaryOpNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst = arg1->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  NodeID arg2subst = arg2->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  NodeID arg3subst = arg3->substituteEndoLagGreaterThanTwo(subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

NodeID
TrinaryOpNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  if (maxExoLead() == 0)
    return const_cast<TrinaryOpNode *>(this);
  else
    return createExoLeadAuxiliaryVarForMyself(subst_table, neweqs);
}

NodeID
TrinaryOpNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  NodeID arg1subst = arg1->substituteExoLag(subst_table, neweqs);
  NodeID arg2subst = arg2->substituteExoLag(subst_table, neweqs);
  NodeID arg3subst = arg3->substituteExoLag(subst_table, neweqs);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

NodeID
TrinaryOpNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  NodeID arg1subst = arg1->substituteExpectation(subst_table, neweqs, partial_information_model);
  NodeID arg2subst = arg2->substituteExpectation(subst_table, neweqs, partial_information_model);
  NodeID arg3subst = arg3->substituteExpectation(subst_table, neweqs, partial_information_model);
  return buildSimilarTrinaryOpNode(arg1subst, arg2subst, arg3subst, datatree);
}

bool
TrinaryOpNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
TrinaryOpNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

UnknownFunctionNode::UnknownFunctionNode(DataTree &datatree_arg,
                                         int symb_id_arg,
                                         const vector<NodeID> &arguments_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  arguments(arguments_arg)
{
}

void
UnknownFunctionNode::prepareForDerivation()
{
  cerr << "UnknownFunctionNode::prepareForDerivation: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::computeDerivative(int deriv_id)
{
  cerr << "UnknownFunctionNode::computeDerivative: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::getChainRuleDerivative(int deriv_id, const map<int, NodeID> &recursive_variables)
{
  cerr << "UnknownFunctionNode::getChainRuleDerivative: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                           temporary_terms_type &temporary_terms,
                                           bool is_matlab) const
{
  cerr << "UnknownFunctionNode::computeTemporaryTerms: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

void
UnknownFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                 const temporary_terms_type &temporary_terms) const
{
  output << datatree.symbol_table.getName(symb_id) << "(";
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeOutput(output, output_type, temporary_terms);
    }
  output << ")";
}

void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                           temporary_terms_type &temporary_terms,
                                           map<NodeID, pair<int, int> > &first_occurence,
                                           int Curr_block,
                                           vector< vector<temporary_terms_type> > &v_temporary_terms,
                                           int equation) const
{
  cerr << "UnknownFunctionNode::computeTemporaryTerms: not implemented" << endl;
  exit(EXIT_FAILURE);
}

void
UnknownFunctionNode::collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const
{
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    (*it)->collectVariables(type_arg, result);
}

void
UnknownFunctionNode::collectTemporary_terms(const temporary_terms_type &temporary_terms, temporary_terms_inuse_type &temporary_terms_inuse, int Curr_Block) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnknownFunctionNode *>(this));
  if (it != temporary_terms.end())
    temporary_terms_inuse.insert(idx);
  else
    {
      //arg->collectTemporary_terms(temporary_terms, result);
    }
}

double
UnknownFunctionNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  throw EvalException();
}

void
UnknownFunctionNode::compile(ostream &CompileCode, bool lhs_rhs, const temporary_terms_type &temporary_terms, map_idx_type &map_idx, bool dynamic, bool steady_dynamic) const
{
  cerr << "UnknownFunctionNode::compile: operation impossible!" << endl;
  exit(EXIT_FAILURE);
}

pair<int, NodeID>
UnknownFunctionNode::normalizeEquation(int var_endo, vector<pair<int, pair<NodeID, NodeID> > >  &List_of_Op_RHS) const
{
  vector<pair<bool, NodeID> > V_arguments;
  vector<NodeID> V_NodeID;
  bool present = false;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    {
      V_arguments.push_back((*it)->normalizeEquation(var_endo, List_of_Op_RHS));
      present = present || V_arguments[V_arguments.size()-1].first;
      V_NodeID.push_back(V_arguments[V_arguments.size()-1].second);
    }
  if (!present)
    return (make_pair(0, datatree.AddUnknownFunction(datatree.symbol_table.getName(symb_id), V_NodeID)));
  else
    return (make_pair(1, (NodeID) NULL));
}

NodeID
UnknownFunctionNode::toStatic(DataTree &static_datatree) const
{
  vector<NodeID> static_arguments;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    static_arguments.push_back((*it)->toStatic(static_datatree));
  return static_datatree.AddUnknownFunction(datatree.symbol_table.getName(symb_id), static_arguments);
}

int
UnknownFunctionNode::maxEndoLead() const
{
  int val = 0;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxEndoLead());
  return val;
}

int
UnknownFunctionNode::maxExoLead() const
{
  int val = 0;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxExoLead());
  return val;
}

int
UnknownFunctionNode::maxEndoLag() const
{
  int val = 0;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxEndoLag());
  return val;
}

int
UnknownFunctionNode::maxExoLag() const
{
  int val = 0;
  for (vector<NodeID>::const_iterator it = arguments.begin();
       it != arguments.end(); it++)
    val = max(val, (*it)->maxExoLag());
  return val;
}

NodeID
UnknownFunctionNode::decreaseLeadsLags(int n) const
{
  cerr << "UnknownFunctionNode::decreaseLeadsLags: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::decreaseLeadsLagsPredeterminedVariables() const
{
  cerr << "UnknownFunctionNode::decreaseLeadsLagsPredeterminedVariables: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "UnknownFunctionNode::substituteEndoLeadGreaterThanTwo: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "UnknownFunctionNode::substituteEndoLagGreaterThanTwo: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "UnknownFunctionNode::substituteExoLead: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const
{
  cerr << "UnknownFunctionNode::substituteExoLag: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

NodeID
UnknownFunctionNode::substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const
{
  cerr << "UnknownFunctionNode::substituteExpectation: not implemented!" << endl;
  exit(EXIT_FAILURE);
}

bool
UnknownFunctionNode::isNumConstNodeEqualTo(double value) const
{
  return false;
}

bool
UnknownFunctionNode::isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const
{
  return false;
}

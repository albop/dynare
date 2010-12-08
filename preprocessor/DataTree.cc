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

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "DataTree.hh"

DataTree::DataTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg) :
  symbol_table(symbol_table_arg),
  num_constants(num_constants_arg),
  node_counter(0)
{
  Zero = AddNumConstant("0");
  One = AddNumConstant("1");
  Two = AddNumConstant("2");

  MinusOne = AddUMinus(One);

  NaN = AddNumConstant("NaN");
  Infinity = AddNumConstant("Inf");
  MinusInfinity = AddUMinus(Infinity);

  Pi = AddNumConstant("3.141592653589793");
}

DataTree::~DataTree()
{
  for (node_list_type::iterator it = node_list.begin(); it != node_list.end(); it++)
    delete *it;
}

NodeID
DataTree::AddNumConstant(const string &value)
{
  int id = num_constants.AddConstant(value);

  num_const_node_map_type::iterator it = num_const_node_map.find(id);
  if (it != num_const_node_map.end())
    return it->second;
  else
    return new NumConstNode(*this, id);
}

VariableNode *
DataTree::AddVariableInternal(int symb_id, int lag)
{
  variable_node_map_type::iterator it = variable_node_map.find(make_pair(symb_id, lag));
  if (it != variable_node_map.end())
    return it->second;
  else
    return new VariableNode(*this, symb_id, lag);
}

VariableNode *
DataTree::AddVariable(const string &name, int lag)
{
  int symb_id = symbol_table.getID(name);
  return AddVariable(symb_id, lag);
}

VariableNode *
DataTree::AddVariable(int symb_id, int lag)
{
  assert(lag == 0);
  return AddVariableInternal(symb_id, lag);
}

NodeID
DataTree::AddPlus(NodeID iArg1, NodeID iArg2)
{
  if (iArg1 != Zero && iArg2 != Zero)
    {
      // Simplify x+(-y) in x-y
      UnaryOpNode *uarg2 = dynamic_cast<UnaryOpNode *>(iArg2);
      if (uarg2 != NULL && uarg2->get_op_code() == oUminus)
        return AddMinus(iArg1, uarg2->get_arg());

      // To treat commutativity of "+"
      // Nodes iArg1 and iArg2 are sorted by index
      if (iArg1->idx > iArg2->idx)
        {
          NodeID tmp = iArg1;
          iArg1 = iArg2;
          iArg2 = tmp;
        }
      return AddBinaryOp(iArg1, oPlus, iArg2);
    }
  else if (iArg1 != Zero)
    return iArg1;
  else if (iArg2 != Zero)
    return iArg2;
  else
    return Zero;
}

NodeID
DataTree::AddMinus(NodeID iArg1, NodeID iArg2)
{
  if (iArg2 == Zero)
    return iArg1;

  if (iArg1 == Zero)
    return AddUMinus(iArg2);

  if (iArg1 == iArg2)
    return Zero;

  return AddBinaryOp(iArg1, oMinus, iArg2);
}

NodeID
DataTree::AddUMinus(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      // Simplify -(-x) in x
      UnaryOpNode *uarg = dynamic_cast<UnaryOpNode *>(iArg1);
      if (uarg != NULL && uarg->get_op_code() == oUminus)
        return uarg->get_arg();

      return AddUnaryOp(oUminus, iArg1);
    }
  else
    return Zero;
}

NodeID
DataTree::AddTimes(NodeID iArg1, NodeID iArg2)
{
  if (iArg1 == MinusOne)
    return AddUMinus(iArg2);
  else if (iArg2 == MinusOne)
    return AddUMinus(iArg1);
  else if (iArg1 != Zero && iArg1 != One && iArg2 != Zero && iArg2 != One)
    {
      // To treat commutativity of "*"
      // Nodes iArg1 and iArg2 are sorted by index
      if (iArg1->idx > iArg2->idx)
        {
          NodeID tmp = iArg1;
          iArg1 = iArg2;
          iArg2 = tmp;
        }
      return AddBinaryOp(iArg1, oTimes, iArg2);
    }
  else if (iArg1 != Zero && iArg1 != One && iArg2 == One)
    return iArg1;
  else if (iArg2 != Zero && iArg2 != One && iArg1 == One)
    return iArg2;
  else if (iArg2 == One && iArg1 == One)
    return One;
  else
    return Zero;
}

NodeID
DataTree::AddDivide(NodeID iArg1, NodeID iArg2)
{
  if (iArg2 == One)
    return iArg1;

  // This test should be before the next two, otherwise 0/0 won't be rejected
  if (iArg2 == Zero)
    {
      cerr << "ERROR: Division by zero!" << endl;
      exit(EXIT_FAILURE);
    }

  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == iArg2)
    return One;

  return AddBinaryOp(iArg1, oDivide, iArg2);
}

NodeID
DataTree::AddLess(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oLess, iArg2);
}

NodeID
DataTree::AddGreater(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oGreater, iArg2);
}

NodeID
DataTree::AddLessEqual(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oLessEqual, iArg2);
}

NodeID
DataTree::AddGreaterEqual(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oGreaterEqual, iArg2);
}

NodeID
DataTree::AddEqualEqual(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oEqualEqual, iArg2);
}

NodeID
DataTree::AddDifferent(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oDifferent, iArg2);
}

NodeID
DataTree::AddPower(NodeID iArg1, NodeID iArg2)
{
  if (iArg1 != Zero && iArg2 != Zero && iArg2 != One)
    return AddBinaryOp(iArg1, oPower, iArg2);
  else if (iArg2 == One)
    return iArg1;
  else if (iArg2 == Zero)
    return One;
  else
    return Zero;
}

NodeID
DataTree::AddExp(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oExp, iArg1);
  else
    return One;
}

NodeID
DataTree::AddLog(NodeID iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    return AddUnaryOp(oLog, iArg1);
  else if (iArg1 == One)
    return Zero;
  else
    {
      cerr << "ERROR: log(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }
}

NodeID
DataTree::AddLog10(NodeID iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    return AddUnaryOp(oLog10, iArg1);
  else if (iArg1 == One)
    return Zero;
  else
    {
      cerr << "ERROR: log10(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }
}

NodeID
DataTree::AddCos(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oCos, iArg1);
  else
    return One;
}

NodeID
DataTree::AddSin(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSin, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddTan(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oTan, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAcos(NodeID iArg1)
{
  if (iArg1 != One)
    return AddUnaryOp(oAcos, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAsin(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAsin, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAtan(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAtan, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddCosh(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oCosh, iArg1);
  else
    return One;
}

NodeID
DataTree::AddSinh(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSinh, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddTanh(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oTanh, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAcosh(NodeID iArg1)
{
  if (iArg1 != One)
    return AddUnaryOp(oAcosh, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAsinh(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAsinh, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddAtanh(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAtanh, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddSqrt(NodeID iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSqrt, iArg1);
  else
    return Zero;
}

NodeID
DataTree::AddMax(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oMax, iArg2);
}

NodeID
DataTree::AddMin(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oMin, iArg2);
}

NodeID
DataTree::AddNormcdf(NodeID iArg1, NodeID iArg2, NodeID iArg3)
{
  return AddTrinaryOp(iArg1, oNormcdf, iArg2, iArg3);
}

NodeID
DataTree::AddSteadyState(NodeID iArg1)
{
  return AddUnaryOp(oSteadyState, iArg1);
}

NodeID
DataTree::AddExpectation(int iArg1, NodeID iArg2)
{
  return AddUnaryOp(oExpectation, iArg2, iArg1);
}

NodeID
DataTree::AddEqual(NodeID iArg1, NodeID iArg2)
{
  return AddBinaryOp(iArg1, oEqual, iArg2);
}

void
DataTree::AddLocalVariable(const string &name, NodeID value) throw (LocalVariableException)
{
  int id = symbol_table.getID(name);

  assert(symbol_table.getType(id) == eModelLocalVariable);

  // Throw an exception if symbol already declared
  map<int, NodeID>::iterator it = local_variables_table.find(id);
  if (it != local_variables_table.end())
    throw LocalVariableException(name);

  local_variables_table[id] = value;
}

NodeID
DataTree::AddUnknownFunction(const string &function_name, const vector<NodeID> &arguments)
{
  int id = symbol_table.getID(function_name);

  assert(symbol_table.getType(id) == eUnknownFunction);

  return new UnknownFunctionNode(*this, id, arguments);
}

bool
DataTree::isSymbolUsed(int symb_id) const
{
  for (variable_node_map_type::const_iterator it = variable_node_map.begin();
       it != variable_node_map.end(); it++)
    if (it->first.first == symb_id)
      return true;

  if (local_variables_table.find(symb_id) != local_variables_table.end())
    return true;

  return false;
}

int
DataTree::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

int
DataTree::getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

bool
DataTree::isUnaryOpUsed(UnaryOpcode opcode) const
{
  for (unary_op_node_map_type::const_iterator it = unary_op_node_map.begin();
       it != unary_op_node_map.end(); it++)
    if (it->first.first.second == opcode)
      return true;

  return false;
}

bool
DataTree::isBinaryOpUsed(BinaryOpcode opcode) const
{
  for (binary_op_node_map_type::const_iterator it = binary_op_node_map.begin();
       it != binary_op_node_map.end(); it++)
    if (it->first.second == opcode)
      return true;

  return false;
}

bool
DataTree::isTrinaryOpUsed(TrinaryOpcode opcode) const
{
  for (trinary_op_node_map_type::const_iterator it = trinary_op_node_map.begin();
       it != trinary_op_node_map.end(); it++)
    if (it->first.second == opcode)
      return true;

  return false;
}

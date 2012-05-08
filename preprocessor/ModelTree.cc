/*
 * Copyright (C) 2003-2011 Dynare Team
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
#include <cmath>
#include <iostream>
#include <fstream>

#include "ModelTree.hh"
#include "MinimumFeedbackSet.hh"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

using namespace boost;
using namespace MFS;

bool
ModelTree::computeNormalization(const jacob_map_t &contemporaneous_jacobian, bool verbose)
{
  const int n = equation_number();

  assert(n == symbol_table.endo_nbr());

  typedef adjacency_list<vecS, vecS, undirectedS> BipartiteGraph;

  /*
    Vertices 0 to n-1 are for endogenous (using type specific ID)
    Vertices n to 2*n-1 are for equations (using equation no.)
  */
  BipartiteGraph g(2 * n);

  // Fill in the graph
  set<pair<int, int> > endo;

  for (jacob_map_t::const_iterator it = contemporaneous_jacobian.begin(); it != contemporaneous_jacobian.end(); it++)
    add_edge(it->first.first + n, it->first.second, g);

  // Compute maximum cardinality matching
  vector<int> mate_map(2*n);

#if 1
  bool check = checked_edmonds_maximum_cardinality_matching(g, &mate_map[0]);
#else // Alternative way to compute normalization, by giving an initial matching using natural normalizations
  fill(mate_map.begin(), mate_map.end(), graph_traits<BipartiteGraph>::null_vertex());

  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      int j = natural_endo2eqs.find(i)->second;

      put(&mate_map[0], i, n+j);
      put(&mate_map[0], n+j, i);
    }

  edmonds_augmenting_path_finder<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type> augmentor(g, &mate_map[0], get(vertex_index, g));
  bool not_maximum_yet = true;
  while (not_maximum_yet)
    {
      not_maximum_yet = augmentor.augment_matching();
    }
  augmentor.get_current_matching(&mate_map[0]);

  bool check = maximum_cardinality_matching_verifier<BipartiteGraph, size_t *, property_map<BipartiteGraph, vertex_index_t>::type>::verify_matching(g, &mate_map[0], get(vertex_index, g));
#endif

  assert(check);

#ifdef DEBUG
  for (int i = 0; i < n; i++)
    cout << "Endogenous " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
         << " matched with equation " << (mate_map[i]-n+1) << endl;
#endif

  // Create the resulting map, by copying the n first elements of mate_map, and substracting n to them
  endo2eq.resize(equation_number());
  transform(mate_map.begin(), mate_map.begin() + n, endo2eq.begin(), bind2nd(minus<int>(), n));

#ifdef DEBUG
  multimap<int, int> natural_endo2eqs;
  computeNormalizedEquations(natural_endo2eqs);

  int n1 = 0, n2 = 0;

  for (int i = 0; i < symbol_table.endo_nbr(); i++)
    {
      if (natural_endo2eqs.count(i) == 0)
        continue;

      n1++;

      pair<multimap<int, int>::const_iterator, multimap<int, int>::const_iterator> x = natural_endo2eqs.equal_range(i);
      if (find_if(x.first, x.second, compose1(bind2nd(equal_to<int>(), endo2eq[i]), select2nd<multimap<int, int>::value_type>())) == x.second)
        cout << "Natural normalization of variable " << symbol_table.getName(symbol_table.getID(eEndogenous, i))
             << " not used." << endl;
      else
        n2++;
    }

  cout << "Used " << n2 << " natural normalizations out of " << n1 << ", for a total of " << n << " equations." << endl;
#endif

  // Check if all variables are normalized
  vector<int>::const_iterator it = find(mate_map.begin(), mate_map.begin() + n, graph_traits<BipartiteGraph>::null_vertex());
  if (it != mate_map.begin() + n)
    {
      if (verbose)
        cerr << "ERROR: Could not normalize the model. Variable "
             << symbol_table.getName(symbol_table.getID(eEndogenous, it - mate_map.begin()))
             << " is not in the maximum cardinality matching." << endl;
      check = false;
    }
  return check;
}

void
ModelTree::computeNonSingularNormalization(jacob_map_t &contemporaneous_jacobian, double cutoff, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian)
{
  bool check = false;

  cout << "Normalizing the model..." << endl;

  int n = equation_number();

  // compute the maximum value of each row of the contemporaneous Jacobian matrix
  //jacob_map normalized_contemporaneous_jacobian;
  jacob_map_t normalized_contemporaneous_jacobian(contemporaneous_jacobian);
  vector<double> max_val(n, 0.0);
  for (jacob_map_t::const_iterator iter = contemporaneous_jacobian.begin(); iter != contemporaneous_jacobian.end(); iter++)
    if (fabs(iter->second) > max_val[iter->first.first])
      max_val[iter->first.first] = fabs(iter->second);

  for (jacob_map_t::iterator iter = normalized_contemporaneous_jacobian.begin(); iter != normalized_contemporaneous_jacobian.end(); iter++)
    iter->second /= max_val[iter->first.first];

  //We start with the highest value of the cutoff and try to normalize the model
  double current_cutoff = 0.99999999;

  int suppressed = 0;
  while (!check && current_cutoff > 1e-19)
    {
      jacob_map_t tmp_normalized_contemporaneous_jacobian;
      int suppress = 0;
      for (jacob_map_t::iterator iter = normalized_contemporaneous_jacobian.begin(); iter != normalized_contemporaneous_jacobian.end(); iter++)
        if (fabs(iter->second) > max(current_cutoff, cutoff))
          tmp_normalized_contemporaneous_jacobian[make_pair(iter->first.first, iter->first.second)] = iter->second;
        else
          suppress++;

      if (suppress != suppressed)
        check = computeNormalization(tmp_normalized_contemporaneous_jacobian, false);
      suppressed = suppress;
      if (!check)
        {
          current_cutoff /= 2;
          // In this last case try to normalize with the complete jacobian
          if (current_cutoff <= 1e-19)
            check = computeNormalization(normalized_contemporaneous_jacobian, false);
        }
    }

  if (!check)
    {
      cout << "Normalization failed with cutoff, trying symbolic normalization..." << endl;
      //if no non-singular normalization can be found, try to find a normalization even with a potential singularity
      jacob_map_t tmp_normalized_contemporaneous_jacobian;
      set<pair<int, int> > endo;
      for (int i = 0; i < n; i++)
        {
          endo.clear();
          equations[i]->collectEndogenous(endo);
          for (set<pair<int, int> >::const_iterator it = endo.begin(); it != endo.end(); it++)
            tmp_normalized_contemporaneous_jacobian[make_pair(i, it->first)] = 1;
        }
      check = computeNormalization(tmp_normalized_contemporaneous_jacobian, true);
      if (check)
        {
          // Update the jacobian matrix
          for (jacob_map_t::const_iterator it = tmp_normalized_contemporaneous_jacobian.begin(); it != tmp_normalized_contemporaneous_jacobian.end(); it++)
            {
              if (static_jacobian.find(make_pair(it->first.first, it->first.second)) == static_jacobian.end())
                static_jacobian[make_pair(it->first.first, it->first.second)] = 0;
              if (dynamic_jacobian.find(make_pair(0, make_pair(it->first.first, it->first.second))) == dynamic_jacobian.end())
                dynamic_jacobian[make_pair(0, make_pair(it->first.first, it->first.second))] = 0;
              if (contemporaneous_jacobian.find(make_pair(it->first.first, it->first.second)) == contemporaneous_jacobian.end())
                contemporaneous_jacobian[make_pair(it->first.first, it->first.second)] = 0;
              if (first_derivatives.find(make_pair(it->first.first, getDerivID(symbol_table.getID(eEndogenous, it->first.second), 0))) == first_derivatives.end())
                first_derivatives[make_pair(it->first.first, getDerivID(symbol_table.getID(eEndogenous, it->first.second), 0))] = Zero;
            }
        }
    }

  if (!check)
    {
      cerr << "No normalization could be computed. Aborting." << endl;
      exit(EXIT_FAILURE);
    }
}

void
ModelTree::computeNormalizedEquations(multimap<int, int> &endo2eqs) const
{
  for (int i = 0; i < equation_number(); i++)
    {
      VariableNode *lhs = dynamic_cast<VariableNode *>(equations[i]->get_arg1());
      if (lhs == NULL)
        continue;

      int symb_id = lhs->get_symb_id();
      if (symbol_table.getType(symb_id) != eEndogenous)
        continue;

      set<pair<int, int> > endo;
      equations[i]->get_arg2()->collectEndogenous(endo);
      if (endo.find(make_pair(symbol_table.getTypeSpecificID(symb_id), 0)) != endo.end())
        continue;

      endo2eqs.insert(make_pair(symbol_table.getTypeSpecificID(symb_id), i));
      cout << "Endogenous " << symbol_table.getName(symb_id) << " normalized in equation " << (i+1) << endl;
    }
}

void
ModelTree::evaluateAndReduceJacobian(const eval_context_t &eval_context, jacob_map_t &contemporaneous_jacobian, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian, double cutoff, bool verbose)
{
  int nb_elements_contemparenous_Jacobian = 0;
  set<pair<int, int> > jacobian_elements_to_delete;
  for (first_derivatives_t::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int deriv_id = it->first.second;
      if (getTypeByDerivID(deriv_id) == eEndogenous)
        {
          expr_t Id = it->second;
          int eq = it->first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          double val = 0;
          try
            {
              val = Id->eval(eval_context);
            }
          catch (ExprNode::EvalExternalFunctionException &e)
            {
              val = 1;
            }
          catch (ExprNode::EvalException &e)
            {
              cerr << "ERROR: evaluation of Jacobian failed for equation " << eq+1 << " and variable " << symbol_table.getName(symb) << "(" << lag << ") [" << symb << "] !" << endl;
              Id->writeOutput(cerr, oMatlabDynamicModelSparse, temporary_terms);
              cerr << endl;
              exit(EXIT_FAILURE);
            }
          if (fabs(val) < cutoff)
            {
              if (verbose)
                cout << "the coefficient related to variable " << var << " with lag " << lag << " in equation " << eq << " is equal to " << val << " and is set to 0 in the incidence matrix (size=" << symbol_table.endo_nbr() << ")" << endl;
              jacobian_elements_to_delete.insert(make_pair(eq, deriv_id));
            }
          else
            {
              if (lag == 0)
                {
                  nb_elements_contemparenous_Jacobian++;
                  contemporaneous_jacobian[make_pair(eq, var)] = val;
                }
              if (static_jacobian.find(make_pair(eq, var)) != static_jacobian.end())
                static_jacobian[make_pair(eq, var)] += val;
              else
                static_jacobian[make_pair(eq, var)] = val;
              dynamic_jacobian[make_pair(lag, make_pair(eq, var))] = Id;
            }
        }
    }

  // Get rid of the elements of the Jacobian matrix below the cutoff
  for (set<pair<int, int> >::const_iterator it = jacobian_elements_to_delete.begin(); it != jacobian_elements_to_delete.end(); it++)
    first_derivatives.erase(*it);

  if (jacobian_elements_to_delete.size() > 0)
    {
      cout << jacobian_elements_to_delete.size() << " elements among " << first_derivatives.size() << " in the incidence matrices are below the cutoff (" << cutoff << ") and are discarded" << endl
           << "The contemporaneous incidence matrix has " << nb_elements_contemparenous_Jacobian << " elements" << endl;
    }
}

void
ModelTree::computePrologueAndEpilogue(const jacob_map_t &static_jacobian_arg, vector<int> &equation_reordered, vector<int> &variable_reordered)
{
  vector<int> eq2endo(equation_number(), 0);
  equation_reordered.resize(equation_number());
  variable_reordered.resize(equation_number());
  bool *IM;
  int n = equation_number();
  IM = (bool *) calloc(n*n, sizeof(bool));
  int i = 0;
  for (vector<int>::const_iterator it = endo2eq.begin(); it != endo2eq.end(); it++, i++)
    {
      eq2endo[*it] = i;
      equation_reordered[i] = i;
      variable_reordered[*it] = i;
    }
  for (jacob_map_t::const_iterator it = static_jacobian_arg.begin(); it != static_jacobian_arg.end(); it++)
    IM[it->first.first * n + endo2eq[it->first.second]] = true;
  bool something_has_been_done = true;
  prologue = 0;
  int k = 0;
  // Find the prologue equations and place first the AR(1) shock equations first
  while (something_has_been_done)
    {
      int tmp_prologue = prologue;
      something_has_been_done = false;
      for (int i = prologue; i < n; i++)
        {
          int nze = 0;
          for (int j = tmp_prologue; j < n; j++)
            if (IM[i * n + j])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[tmp_prologue * n + j];
                  IM[tmp_prologue * n + j] = IM[i * n + j];
                  IM[i * n + j] = tmp_bool;
                }
              int tmp = equation_reordered[tmp_prologue];
              equation_reordered[tmp_prologue] = equation_reordered[i];
              equation_reordered[i] = tmp;
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[j * n + tmp_prologue];
                  IM[j * n + tmp_prologue] = IM[j * n + k];
                  IM[j * n + k] = tmp_bool;
                }
              tmp = variable_reordered[tmp_prologue];
              variable_reordered[tmp_prologue] = variable_reordered[k];
              variable_reordered[k] = tmp;
              tmp_prologue++;
              something_has_been_done = true;
            }
        }
      prologue = tmp_prologue;
    }

  something_has_been_done = true;
  epilogue = 0;
  // Find the epilogue equations
  while (something_has_been_done)
    {
      int tmp_epilogue = epilogue;
      something_has_been_done = false;
      for (int i = prologue; i < n - (int) epilogue; i++)
        {
          int nze = 0;
          for (int j = prologue; j < n - tmp_epilogue; j++)
            if (IM[j * n + i])
              {
                nze++;
                k = j;
              }
          if (nze == 1)
            {
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[(n - 1 - tmp_epilogue) * n + j];
                  IM[(n - 1 - tmp_epilogue) * n + j] = IM[k * n + j];
                  IM[k * n + j] = tmp_bool;
                }
              int tmp = equation_reordered[n - 1 - tmp_epilogue];
              equation_reordered[n - 1 - tmp_epilogue] = equation_reordered[k];
              equation_reordered[k] = tmp;
              for (int j = 0; j < n; j++)
                {
                  bool tmp_bool = IM[j * n + n - 1 - tmp_epilogue];
                  IM[j * n + n - 1 - tmp_epilogue] = IM[j * n + i];
                  IM[j * n + i] = tmp_bool;
                }
              tmp = variable_reordered[n - 1 - tmp_epilogue];
              variable_reordered[n - 1 - tmp_epilogue] = variable_reordered[i];
              variable_reordered[i] = tmp;
              tmp_epilogue++;
              something_has_been_done = true;
            }
        }
      epilogue = tmp_epilogue;
    }
  free(IM);
}

equation_type_and_normalized_equation_t
ModelTree::equationTypeDetermination(const map<pair<int, pair<int, int> >, expr_t> &first_order_endo_derivatives, const vector<int> &Index_Var_IM, const vector<int> &Index_Equ_IM, int mfs) const
{
  expr_t lhs, rhs;
  BinaryOpNode *eq_node;
  EquationType Equation_Simulation_Type;
  equation_type_and_normalized_equation_t V_Equation_Simulation_Type(equations.size());
  for (unsigned int i = 0; i < equations.size(); i++)
    {
      int eq = Index_Equ_IM[i];
      int var = Index_Var_IM[i];
      eq_node = equations[eq];
      lhs = eq_node->get_arg1();
      rhs = eq_node->get_arg2();
      Equation_Simulation_Type = E_SOLVE;
      map<pair<int, pair<int, int> >, expr_t>::const_iterator derivative = first_order_endo_derivatives.find(make_pair(eq, make_pair(var, 0)));
      pair<bool, expr_t> res;
      if (derivative != first_order_endo_derivatives.end())
        {
          set<pair<int, int> > result;
          derivative->second->collectEndogenous(result);
          set<pair<int, int> >::const_iterator d_endo_variable = result.find(make_pair(var, 0));
          //Determine whether the equation could be evaluated rather than to be solved
          if (lhs->isVariableNodeEqualTo(eEndogenous, Index_Var_IM[i], 0) && derivative->second->isNumConstNodeEqualTo(1))
            {
              Equation_Simulation_Type = E_EVALUATE;
            }
          else
            {
              vector<pair<int, pair<expr_t, expr_t> > > List_of_Op_RHS;
              res =  equations[eq]->normalizeEquation(var, List_of_Op_RHS);
              if (mfs == 2)
                {
                  if (d_endo_variable == result.end() && res.second)
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
              else if (mfs == 3)
                {
                  if (res.second) // The equation could be solved analytically
                    Equation_Simulation_Type = E_EVALUATE_S;
                }
            }
        }
      V_Equation_Simulation_Type[eq] = make_pair(Equation_Simulation_Type, dynamic_cast<BinaryOpNode *>(res.second));
    }
  return (V_Equation_Simulation_Type);
}

void
ModelTree::getVariableLeadLagByBlock(const dynamic_jacob_map_t &dynamic_jacobian, const vector<int> &components_set, int nb_blck_sim, lag_lead_vector_t &equation_lead_lag, lag_lead_vector_t &variable_lead_lag, const vector<int> &equation_reordered, const vector<int> &variable_reordered) const
{
  int nb_endo = symbol_table.endo_nbr();
  variable_lead_lag = lag_lead_vector_t(nb_endo, make_pair(0, 0));
  equation_lead_lag = lag_lead_vector_t(nb_endo, make_pair(0, 0));
  vector<int> variable_blck(nb_endo), equation_blck(nb_endo);
  for (int i = 0; i < nb_endo; i++)
    {
      if (i < (int) prologue)
        {
          variable_blck[variable_reordered[i]] = i;
          equation_blck[equation_reordered[i]] = i;
        }
      else if (i < (int) (components_set.size() + prologue))
        {
          variable_blck[variable_reordered[i]] = components_set[i-prologue] + prologue;
          equation_blck[equation_reordered[i]] = components_set[i-prologue] + prologue;
        }
      else
        {
          variable_blck[variable_reordered[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
          equation_blck[equation_reordered[i]] = i- (nb_endo - nb_blck_sim - prologue - epilogue);
        }
    }
  for (dynamic_jacob_map_t::const_iterator it = dynamic_jacobian.begin(); it != dynamic_jacobian.end(); it++)
    {
      int lag = it->first.first;
      int j_1 = it->first.second.first;
      int i_1 = it->first.second.second;
      if (variable_blck[i_1] == equation_blck[j_1])
        {
          if (lag > variable_lead_lag[i_1].second)
            variable_lead_lag[i_1] = make_pair(variable_lead_lag[i_1].first, lag);
          if (lag < -variable_lead_lag[i_1].first)
            variable_lead_lag[i_1] = make_pair(-lag, variable_lead_lag[i_1].second);
          if (lag > equation_lead_lag[j_1].second)
            equation_lead_lag[j_1] = make_pair(equation_lead_lag[j_1].first, lag);
          if (lag < -equation_lead_lag[j_1].first)
            equation_lead_lag[j_1] = make_pair(-lag, equation_lead_lag[j_1].second);
        }
    }
}

void
ModelTree::computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered, vector<pair<int, int> > &blocks, const equation_type_and_normalized_equation_t &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered, lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed) const
{
  int nb_var = variable_reordered.size();
  int n = nb_var - prologue - epilogue;

  AdjacencyList_t G2(n);

  // It is necessary to manually initialize vertex_index property since this graph uses listS and not vecS as underlying vertex container
  property_map<AdjacencyList_t, vertex_index_t>::type v_index = get(vertex_index, G2);
  for (int i = 0; i < n; i++)
    put(v_index, vertex(i, G2), i);

  vector<int> reverse_equation_reordered(nb_var), reverse_variable_reordered(nb_var);

  for (int i = 0; i < nb_var; i++)
    {
      reverse_equation_reordered[equation_reordered[i]] = i;
      reverse_variable_reordered[variable_reordered[i]] = i;
    }

  for (jacob_map_t::const_iterator it = static_jacobian.begin(); it != static_jacobian.end(); it++)
    if (reverse_equation_reordered[it->first.first] >= (int) prologue && reverse_equation_reordered[it->first.first] < (int) (nb_var - epilogue)
        && reverse_variable_reordered[it->first.second] >= (int) prologue && reverse_variable_reordered[it->first.second] < (int) (nb_var - epilogue)
        && it->first.first != endo2eq[it->first.second])
      add_edge(vertex(reverse_equation_reordered[endo2eq[it->first.second]]-prologue, G2),
               vertex(reverse_equation_reordered[it->first.first]-prologue, G2),
               G2);

  vector<int> endo2block(num_vertices(G2)), discover_time(num_vertices(G2));
  iterator_property_map<int *, property_map<AdjacencyList_t, vertex_index_t>::type, int, int &> endo2block_map(&endo2block[0], get(vertex_index, G2));

  // Compute strongly connected components
  int num = strong_components(G2, endo2block_map);

  blocks = vector<pair<int, int> >(num, make_pair(0, 0));

  // Create directed acyclic graph associated to the strongly connected components
  typedef adjacency_list<vecS, vecS, directedS> DirectedGraph;
  DirectedGraph dag(num);

  for (unsigned int i = 0; i < num_vertices(G2); i++)
    {
      AdjacencyList_t::out_edge_iterator it_out, out_end;
      AdjacencyList_t::vertex_descriptor vi = vertex(i, G2);
      for (tie(it_out, out_end) = out_edges(vi, G2); it_out != out_end; ++it_out)
        {
          int t_b = endo2block_map[target(*it_out, G2)];
          int s_b = endo2block_map[source(*it_out, G2)];
          if (s_b != t_b)
            add_edge(s_b, t_b, dag);
        }
    }

  // Compute topological sort of DAG (ordered list of unordered SCC)
  deque<int> ordered2unordered;
  topological_sort(dag, front_inserter(ordered2unordered)); // We use a front inserter because topological_sort returns the inverse order

  // Construct mapping from unordered SCC to ordered SCC
  vector<int> unordered2ordered(num);
  for (int i = 0; i < num; i++)
    unordered2ordered[ordered2unordered[i]] = i;

  //This vector contains for each block:
  //   - first set = equations belonging to the block,
  //   - second set = the feeback variables,
  //   - third vector = the reordered non-feedback variables.
  vector<pair<set<int>, pair<set<int>, vector<int> > > > components_set(num);
  for (unsigned int i = 0; i < endo2block.size(); i++)
    {
      endo2block[i] = unordered2ordered[endo2block[i]];
      blocks[endo2block[i]].first++;
      components_set[endo2block[i]].first.insert(i);
    }

  getVariableLeadLagByBlock(dynamic_jacobian, endo2block, num, equation_lag_lead, variable_lag_lead, equation_reordered, variable_reordered);

  vector<int> tmp_equation_reordered(equation_reordered), tmp_variable_reordered(variable_reordered);
  int order = prologue;
  //Add a loop on vertices which could not be normalized or vertices related to lead variables => force those vertices to belong to the feedback set
  if (select_feedback_variable)
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE
            || variable_lag_lead[variable_reordered[i+prologue]].second > 0
            || variable_lag_lead[variable_reordered[i+prologue]].first > 0
            || equation_lag_lead[equation_reordered[i+prologue]].second > 0
            || equation_lag_lead[equation_reordered[i+prologue]].first > 0
            || mfs == 0)
          add_edge(vertex(i, G2), vertex(i, G2), G2);
    }
  else
    {
      for (int i = 0; i < n; i++)
        if (Equation_Type[equation_reordered[i+prologue]].first == E_SOLVE || mfs == 0)
          add_edge(vertex(i, G2), vertex(i, G2), G2);
    }
  //Determines the dynamic structure of each equation
  n_static = vector<unsigned int>(prologue+num+epilogue, 0);
  n_forward = vector<unsigned int>(prologue+num+epilogue, 0);
  n_backward = vector<unsigned int>(prologue+num+epilogue, 0);
  n_mixed = vector<unsigned int>(prologue+num+epilogue, 0);

  for (int i = 0; i < (int) prologue; i++)
    {
      if      (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
        n_mixed[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second != 0)
        n_forward[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first != 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
        n_backward[i]++;
      else if (variable_lag_lead[tmp_variable_reordered[i]].first == 0 && variable_lag_lead[tmp_variable_reordered[i]].second == 0)
        n_static[i]++;
    }
  //For each block, the minimum set of feedback variable is computed
  // and the non-feedback variables are reordered to get
  // a sub-recursive block without feedback variables

  for (int i = 0; i < num; i++)
    {
      AdjacencyList_t G = extract_subgraph(G2, components_set[i].first);
      set<int> feed_back_vertices;
      //Print(G);
      AdjacencyList_t G1 = Minimal_set_of_feedback_vertex(feed_back_vertices, G);
      property_map<AdjacencyList_t, vertex_index_t>::type v_index = get(vertex_index, G);
      components_set[i].second.first = feed_back_vertices;
      blocks[i].second = feed_back_vertices.size();
      vector<int> Reordered_Vertice;
      Reorder_the_recursive_variables(G, feed_back_vertices, Reordered_Vertice);

      //First we have the recursive equations conditional on feedback variables
      for (int j = 0; j < 4; j++)
        {
          for (vector<int>::iterator its = Reordered_Vertice.begin(); its != Reordered_Vertice.end(); its++)
            {
              bool something_done = false;
              if      (j == 2 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].second != 0)
                {
                  n_mixed[prologue+i]++;
                  something_done = true;
                }
              else if (j == 3 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].second != 0)
                {
                  n_forward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 1 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].second == 0)
                {
                  n_backward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 0 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[*its +prologue]].second == 0)
                {
                  n_static[prologue+i]++;
                  something_done = true;
                }
              if (something_done)
                {
                  equation_reordered[order] = tmp_equation_reordered[*its+prologue];
                  variable_reordered[order] = tmp_variable_reordered[*its+prologue];
                  order++;
                }
            }
        }
      components_set[i].second.second = Reordered_Vertice;
      //Second we have the equations related to the feedback variables
      for (int j = 0; j < 4; j++)
        {
          for (set<int>::iterator its = feed_back_vertices.begin(); its != feed_back_vertices.end(); its++)
            {
              bool something_done = false;
              if      (j == 2 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].second != 0)
                {
                  n_mixed[prologue+i]++;
                  something_done = true;
                }
              else if (j == 3 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].second != 0)
                {
                  n_forward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 1 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].first != 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].second == 0)
                {
                  n_backward[prologue+i]++;
                  something_done = true;
                }
              else if (j == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].first == 0 && variable_lag_lead[tmp_variable_reordered[v_index[vertex(*its, G)]+prologue]].second == 0)
                {
                  n_static[prologue+i]++;
                  something_done = true;
                }
              if (something_done)
                {
                  equation_reordered[order] = tmp_equation_reordered[v_index[vertex(*its, G)]+prologue];
                  variable_reordered[order] = tmp_variable_reordered[v_index[vertex(*its, G)]+prologue];
                  order++;
                }
            }
        }
    }

  for (int i = 0; i < (int) epilogue; i++)
    {
      if      (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
        n_mixed[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second != 0)
        n_forward[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first != 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
        n_backward[prologue+num+i]++;
      else if (variable_lag_lead[tmp_variable_reordered[prologue+n+i]].first == 0 && variable_lag_lead[tmp_variable_reordered[prologue+n+i]].second == 0)
        n_static[prologue+num+i]++;
    }

  inv_equation_reordered = vector<int>(nb_var);
  inv_variable_reordered = vector<int>(nb_var);
  for (int i = 0; i < nb_var; i++)
    {
      inv_variable_reordered[variable_reordered[i]] = i;
      inv_equation_reordered[equation_reordered[i]] = i;
    }
}

void
ModelTree::printBlockDecomposition(const vector<pair<int, int> > &blocks) const
{
  int largest_block = 0;
  int Nb_SimulBlocks = 0;
  int Nb_feedback_variable = 0;
  unsigned int Nb_TotalBlocks = getNbBlocks();
  for (unsigned int block = 0; block < Nb_TotalBlocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      if (simulation_type == SOLVE_FORWARD_COMPLETE || simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          Nb_SimulBlocks++;
          int size = getBlockSize(block);
          if (size > largest_block)
            {
              largest_block = size;
              Nb_feedback_variable = blocks[Nb_SimulBlocks-1].second;
            }
        }
    }

  int Nb_RecursBlocks = Nb_TotalBlocks - Nb_SimulBlocks;
  cout << Nb_TotalBlocks << " block(s) found:" << endl
       << "  " << Nb_RecursBlocks << " recursive block(s) and " << Nb_SimulBlocks << " simultaneous block(s)." << endl
       << "  the largest simultaneous block has " << largest_block << " equation(s)" << endl
       << "                                 and " << Nb_feedback_variable << " feedback variable(s)." << endl;
}

block_type_firstequation_size_mfs_t
ModelTree::reduceBlocksAndTypeDetermination(const dynamic_jacob_map_t &dynamic_jacobian, vector<pair<int, int> > &blocks, const equation_type_and_normalized_equation_t &Equation_Type, const vector<int> &variable_reordered, const vector<int> &equation_reordered, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed, vector<pair< pair<int, int>, pair<int, int> > > &block_col_type)
{
  int i = 0;
  int count_equ = 0, blck_count_simult = 0;
  int Blck_Size, MFS_Size;
  int Lead, Lag;
  block_type_firstequation_size_mfs_t block_type_size_mfs;
  BlockSimulationType Simulation_Type, prev_Type = UNKNOWN;
  int eq = 0;
  unsigned int l_n_static = 0;
  unsigned int l_n_forward = 0;
  unsigned int l_n_backward = 0;
  unsigned int l_n_mixed = 0;
  for (i = 0; i < (int) (prologue+blocks.size()+epilogue); i++)
    {
      int first_count_equ = count_equ;
      if (i < (int) prologue)
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }
      else if (i < (int) (prologue+blocks.size()))
        {
          Blck_Size = blocks[blck_count_simult].first;
          MFS_Size = blocks[blck_count_simult].second;
          blck_count_simult++;
        }
      else if (i < (int) (prologue+blocks.size()+epilogue))
        {
          Blck_Size = 1;
          MFS_Size = 1;
        }

      Lag = Lead = 0;
      set<pair<int, int> > endo;
      for (count_equ  = first_count_equ; count_equ  < Blck_Size+first_count_equ; count_equ++)
        {
          endo.clear();
          equations[equation_reordered[count_equ]]->collectEndogenous(endo);
          for (set<pair<int, int> >::const_iterator it = endo.begin(); it != endo.end(); it++)
            {
              int curr_variable = it->first;
              int curr_lag = it->second;
              vector<int>::const_iterator it2 = find(variable_reordered.begin()+first_count_equ, variable_reordered.begin()+(first_count_equ+Blck_Size), curr_variable);
              if (it2 != variable_reordered.begin()+(first_count_equ+Blck_Size))
                if (dynamic_jacobian.find(make_pair(curr_lag, make_pair(equation_reordered[count_equ], curr_variable))) != dynamic_jacobian.end())
                  {
                    if (curr_lag > Lead)
                      Lead = curr_lag;
                    else if (-curr_lag > Lag)
                      Lag = -curr_lag;
                  }
            }
        }
      if ((Lag > 0) && (Lead > 0))
        {
          if (Blck_Size == 1)
            Simulation_Type = SOLVE_TWO_BOUNDARIES_SIMPLE;
          else
            Simulation_Type = SOLVE_TWO_BOUNDARIES_COMPLETE;
        }
      else if (Blck_Size > 1)
        {
          if (Lead > 0)
            Simulation_Type = SOLVE_BACKWARD_COMPLETE;
          else
            Simulation_Type = SOLVE_FORWARD_COMPLETE;
        }
      else
        {
          if (Lead > 0)
            Simulation_Type = SOLVE_BACKWARD_SIMPLE;
          else
            Simulation_Type = SOLVE_FORWARD_SIMPLE;
        }
      l_n_static = n_static[i];
      l_n_forward = n_forward[i];
      l_n_backward = n_backward[i];
      l_n_mixed = n_mixed[i];
      if (Blck_Size == 1)
        {
          if (Equation_Type[equation_reordered[eq]].first == E_EVALUATE || Equation_Type[equation_reordered[eq]].first == E_EVALUATE_S)
            {
              if (Simulation_Type == SOLVE_BACKWARD_SIMPLE)
                Simulation_Type = EVALUATE_BACKWARD;
              else if (Simulation_Type == SOLVE_FORWARD_SIMPLE)
                Simulation_Type = EVALUATE_FORWARD;
            }
          if (i > 0)
            {
              bool is_lead = false, is_lag = false;
              int c_Size = (block_type_size_mfs[block_type_size_mfs.size()-1]).second.first;
              int first_equation = (block_type_size_mfs[block_type_size_mfs.size()-1]).first.second;
              if (c_Size > 0 && ((prev_Type ==  EVALUATE_FORWARD && Simulation_Type == EVALUATE_FORWARD && !is_lead)
                  || (prev_Type ==  EVALUATE_BACKWARD && Simulation_Type == EVALUATE_BACKWARD && !is_lag)))
                {
                  for (int j = first_equation; j < first_equation+c_Size; j++)
                    {
                      dynamic_jacob_map_t::const_iterator it = dynamic_jacobian.find(make_pair(-1, make_pair(equation_reordered[eq], variable_reordered[j])));
                      if (it != dynamic_jacobian.end())
                        is_lag = true;
                      it = dynamic_jacobian.find(make_pair(+1, make_pair(equation_reordered[eq], variable_reordered[j])));
                      if (it != dynamic_jacobian.end())
                        is_lead = true;
                    }
                }
              if ((prev_Type ==  EVALUATE_FORWARD && Simulation_Type == EVALUATE_FORWARD && !is_lead)
                  || (prev_Type ==  EVALUATE_BACKWARD && Simulation_Type == EVALUATE_BACKWARD && !is_lag))
                {
                  //merge the current block with the previous one
                  BlockSimulationType c_Type = (block_type_size_mfs[block_type_size_mfs.size()-1]).first.first;
                  c_Size++;
                  block_type_size_mfs[block_type_size_mfs.size()-1] = make_pair(make_pair(c_Type, first_equation), make_pair(c_Size, c_Size));
                  if (block_lag_lead[block_type_size_mfs.size()-1].first > Lag)
                    Lag = block_lag_lead[block_type_size_mfs.size()-1].first;
                  if (block_lag_lead[block_type_size_mfs.size()-1].second > Lead)
                    Lead = block_lag_lead[block_type_size_mfs.size()-1].second;
                  block_lag_lead[block_type_size_mfs.size()-1] = make_pair(Lag, Lead);
                  pair< pair< unsigned int, unsigned int>, pair<unsigned int, unsigned int> > tmp = block_col_type[block_col_type.size()-1];
                  block_col_type[block_col_type.size()-1] = make_pair(make_pair(tmp.first.first+l_n_static, tmp.first.second+l_n_forward), make_pair(tmp.second.first+l_n_backward, tmp.second.second+l_n_mixed));
                }
              else
                {
                  block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
                  block_lag_lead.push_back(make_pair(Lag, Lead));
                  block_col_type.push_back(make_pair(make_pair(l_n_static, l_n_forward), make_pair(l_n_backward, l_n_mixed)));
                }
            }
          else
            {
              block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
              block_lag_lead.push_back(make_pair(Lag, Lead));
              block_col_type.push_back(make_pair(make_pair(l_n_static, l_n_forward), make_pair(l_n_backward, l_n_mixed)));
            }
        }
      else
        {
          block_type_size_mfs.push_back(make_pair(make_pair(Simulation_Type, eq), make_pair(Blck_Size, MFS_Size)));
          block_lag_lead.push_back(make_pair(Lag, Lead));
          block_col_type.push_back(make_pair(make_pair(l_n_static, l_n_forward), make_pair(l_n_backward, l_n_mixed)));
        }
      prev_Type = Simulation_Type;
      eq += Blck_Size;
    }
  return (block_type_size_mfs);
}

vector<bool>
ModelTree::BlockLinear(const blocks_derivatives_t &blocks_derivatives, const vector<int> &variable_reordered) const
{
  unsigned int nb_blocks = getNbBlocks();
  vector<bool> blocks_linear(nb_blocks, true);
  for (unsigned int block = 0; block < nb_blocks; block++)
    {
      BlockSimulationType simulation_type = getBlockSimulationType(block);
      int block_size = getBlockSize(block);
      block_derivatives_equation_variable_laglead_nodeid_t derivatives_block = blocks_derivatives[block];
      int first_variable_position = getBlockFirstEquation(block);
      if (simulation_type == SOLVE_BACKWARD_COMPLETE || simulation_type == SOLVE_FORWARD_COMPLETE)
        {
          for (block_derivatives_equation_variable_laglead_nodeid_t::const_iterator it = derivatives_block.begin(); it != derivatives_block.end(); it++)
            {
              int lag = it->second.first;
              if (lag == 0)
                {
                  expr_t Id = it->second.second;
                  set<pair<int, int> > endogenous;
                  Id->collectEndogenous(endogenous);
                  if (endogenous.size() > 0)
                    {
                      for (int l = 0; l < block_size; l++)
                        {
                          if (endogenous.find(make_pair(variable_reordered[first_variable_position+l], 0)) != endogenous.end())
                            {
                              blocks_linear[block] = false;
                              goto the_end;
                            }
                        }
                    }
                }
            }
        }
      else if (simulation_type == SOLVE_TWO_BOUNDARIES_COMPLETE || simulation_type == SOLVE_TWO_BOUNDARIES_SIMPLE)
        {
          for (block_derivatives_equation_variable_laglead_nodeid_t::const_iterator it = derivatives_block.begin(); it != derivatives_block.end(); it++)
            {
              int lag = it->second.first;
              expr_t Id = it->second.second; //
              set<pair<int, int> > endogenous;
              Id->collectEndogenous(endogenous);
              if (endogenous.size() > 0)
                {
                  for (int l = 0; l < block_size; l++)
                    {
                      if (endogenous.find(make_pair(variable_reordered[first_variable_position+l], lag)) != endogenous.end())
                        {
                          blocks_linear[block] = false;
                          goto the_end;
                        }
                    }
                }
            }
        }
    the_end:
      ;
    }
  return (blocks_linear);
}

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg,
                     ExternalFunctionsTable &external_functions_table_arg) :
  DataTree(symbol_table_arg, num_constants_arg, external_functions_table_arg)
{
  for (int i = 0; i < 3; i++)
    NNZDerivatives[i] = 0;
}

int
ModelTree::equation_number() const
{
  return (equations.size());
}

void
ModelTree::writeDerivative(ostream &output, int eq, int symb_id, int lag,
                           ExprNodeOutputType output_type,
                           const temporary_terms_t &temporary_terms) const
{
  first_derivatives_t::const_iterator it = first_derivatives.find(make_pair(eq, getDerivID(symb_id, lag)));
  if (it != first_derivatives.end())
    (it->second)->writeOutput(output, output_type, temporary_terms);
  else
    output << 0;
}

void
ModelTree::computeJacobian(const set<int> &vars)
{
  for (set<int>::const_iterator it = vars.begin();
       it != vars.end(); it++)
    for (int eq = 0; eq < (int) equations.size(); eq++)
      {
        expr_t d1 = equations[eq]->getDerivative(*it);
        if (d1 == Zero)
          continue;
        first_derivatives[make_pair(eq, *it)] = d1;
        ++NNZDerivatives[0];
      }
}

void
ModelTree::computeHessian(const set<int> &vars)
{
  for (first_derivatives_t::const_iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var1 = it->first.second;
      expr_t d1 = it->second;

      // Store only second derivatives with var2 <= var1
      for (set<int>::const_iterator it2 = vars.begin();
           it2 != vars.end(); it2++)
        {
          int var2 = *it2;
          if (var2 > var1)
            continue;

          expr_t d2 = d1->getDerivative(var2);
          if (d2 == Zero)
            continue;
          second_derivatives[make_pair(eq, make_pair(var1, var2))] = d2;
          if (var2 == var1)
            ++NNZDerivatives[1];
          else
            NNZDerivatives[1] += 2;
        }
    }
}

void
ModelTree::computeThirdDerivatives(const set<int> &vars)
{
  for (second_derivatives_t::const_iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    {
      int eq = it->first.first;

      int var1 = it->first.second.first;
      int var2 = it->first.second.second;
      // By construction, var2 <= var1

      expr_t d2 = it->second;

      // Store only third derivatives such that var3 <= var2 <= var1
      for (set<int>::const_iterator it2 = vars.begin();
           it2 != vars.end(); it2++)
        {
          int var3 = *it2;
          if (var3 > var2)
            continue;

          expr_t d3 = d2->getDerivative(var3);
          if (d3 == Zero)
            continue;
          third_derivatives[make_pair(eq, make_pair(var1, make_pair(var2, var3)))] = d3;
          if (var3 == var2 && var2 == var1)
            ++NNZDerivatives[2];
          else if (var3 == var2 || var2 == var1)
            NNZDerivatives[2] += 3;
          else
            NNZDerivatives[2] += 6;
        }
    }
}

void
ModelTree::computeTemporaryTerms(bool is_matlab)
{
  map<expr_t, int> reference_count;
  temporary_terms.clear();

  for (vector<BinaryOpNode *>::iterator it = equations.begin();
       it != equations.end(); it++)
    (*it)->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (first_derivatives_t::iterator it = first_derivatives.begin();
       it != first_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (second_derivatives_t::iterator it = second_derivatives.begin();
       it != second_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);

  for (third_derivatives_t::iterator it = third_derivatives.begin();
       it != third_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
}

void
ModelTree::writeTemporaryTerms(const temporary_terms_t &tt, ostream &output,
                               ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_t tt2;

  for (temporary_terms_t::const_iterator it = tt.begin();
       it != tt.end(); it++)
    {
      if (dynamic_cast<ExternalFunctionNode *>(*it) != NULL)
        (*it)->writeExternalFunctionOutput(output, output_type, tt2, tef_terms);

      if (IS_C(output_type))
        output << "double ";

      (*it)->writeOutput(output, output_type, tt, tef_terms);
      output << " = ";
      (*it)->writeOutput(output, output_type, tt2, tef_terms);

      if (IS_C(output_type))
        output << ";" << endl;

      // Insert current node into tt2
      tt2.insert(*it);

      if (IS_MATLAB(output_type))
        output << ";" << endl;
    }
}

void
ModelTree::compileTemporaryTerms(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, map_idx_t map_idx, bool dynamic, bool steady_dynamic) const
{
  // Local var used to keep track of temp nodes already written
  temporary_terms_t tt2;
  // To store the functions that have already been written in the form TEF* = ext_fun();
  deriv_node_temp_terms_t tef_terms;
  for (temporary_terms_t::const_iterator it = tt.begin();
       it != tt.end(); it++)
    {
      if (dynamic_cast<ExternalFunctionNode *>(*it) != NULL)
        {
          (*it)->compileExternalFunctionOutput(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
        }

      FNUMEXPR_ fnumexpr(TemporaryTerm, (int) (map_idx.find((*it)->idx)->second));
      fnumexpr.write(code_file, instruction_number);
      (*it)->compile(code_file, instruction_number, false, tt2, map_idx, dynamic, steady_dynamic, tef_terms);
      if (dynamic)
        {
          FSTPT_ fstpt((int) (map_idx.find((*it)->idx)->second));
          fstpt.write(code_file, instruction_number);
        }
      else
        {
          FSTPST_ fstpst((int) (map_idx.find((*it)->idx)->second));
          fstpst.write(code_file, instruction_number);
        }
      // Insert current node into tt2
      tt2.insert(*it);
    }
}

void
ModelTree::writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const
{
  /* Collect all model local variables appearing in equations, and print only
     them. Printing unused model local variables can lead to a crash (see
     ticket #101). */
  set<int> used_local_vars;

  // Use an empty set for the temporary terms
  const temporary_terms_t tt;

  for (size_t i = 0; i < equations.size(); i++)
    equations[i]->collectModelLocalVariables(used_local_vars);

  for (set<int>::const_iterator it = used_local_vars.begin();
       it != used_local_vars.end(); ++it)
    {
      int id = *it;
      expr_t value = local_variables_table.find(id)->second;
      value->writeExternalFunctionOutput(output, output_type, tt, tef_terms);

      if (IS_C(output_type))
        output << "double ";

      /* We append underscores to avoid name clashes with "g1" or "oo_" (see
         also VariableNode::writeOutput) */
      output << symbol_table.getName(id) << "__ = ";
      value->writeOutput(output, output_type, tt, tef_terms);
      output << ";" << endl;
    }
}

void
ModelTree::writeModelEquations(ostream &output, ExprNodeOutputType output_type) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->get_arg1();
      expr_t rhs = eq_node->get_arg2();

      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_t());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          output << "lhs =";
          lhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;

          output << "rhs =";
          rhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;

          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << "= lhs-rhs;" << endl;
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          output << "residual" << LEFT_ARRAY_SUBSCRIPT(output_type)
                 << eq + ARRAY_SUBSCRIPT_OFFSET(output_type)
                 << RIGHT_ARRAY_SUBSCRIPT(output_type)
                 << " = ";
          lhs->writeOutput(output, output_type, temporary_terms);
          output << ";" << endl;
        }
    }
}

void
ModelTree::compileModelEquations(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const
{
  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];
      expr_t lhs = eq_node->get_arg1();
      expr_t rhs = eq_node->get_arg2();
      FNUMEXPR_ fnumexpr(ModelEquation, eq);
      fnumexpr.write(code_file, instruction_number);
      // Test if the right hand side of the equation is empty.
      double vrhs = 1.0;
      try
        {
          vrhs = rhs->eval(eval_context_t());
        }
      catch (ExprNode::EvalException &e)
        {
        }

      if (vrhs != 0) // The right hand side of the equation is not empty ==> residual=lhs-rhs;
        {
          lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          rhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);

          FBINARY_ fbinary(oMinus);
          fbinary.write(code_file, instruction_number);

          FSTPR_ fstpr(eq);
          fstpr.write(code_file, instruction_number);
        }
      else // The right hand side of the equation is empty ==> residual=lhs;
        {
          lhs->compile(code_file, instruction_number, false, temporary_terms, map_idx, dynamic, steady_dynamic);
          FSTPR_ fstpr(eq);
          fstpr.write(code_file, instruction_number);
        }
    }
}

void
ModelTree::Write_Inf_To_Bin_File(const string &basename,
                                 int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const
{
  int j;
  std::ofstream SaveCode;
  const string bin_basename = basename + ".bin";
  if (file_open)
    SaveCode.open(bin_basename.c_str(), ios::out | ios::in | ios::binary | ios::ate);
  else
    SaveCode.open(bin_basename.c_str(), ios::out | ios::binary);
  if (!SaveCode.is_open())
    {
      cout << "Error : Can't open file \"" << bin_basename << "\" for writing\n";
      exit(EXIT_FAILURE);
    }
  u_count_int = 0;
  for (first_derivatives_t::const_iterator it = first_derivatives.begin(); it != first_derivatives.end(); it++)
    {
      int deriv_id = it->first.second;
      if (getTypeByDerivID(deriv_id) == eEndogenous)
        {
          int eq = it->first.first;
          int symb = getSymbIDByDerivID(deriv_id);
          int var = symbol_table.getTypeSpecificID(symb);
          int lag = getLagByDerivID(deriv_id);
          SaveCode.write(reinterpret_cast<char *>(&eq), sizeof(eq));
          int varr = var + lag * block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&varr), sizeof(varr));
          SaveCode.write(reinterpret_cast<char *>(&lag), sizeof(lag));
          int u = u_count_int + block_mfs;
          SaveCode.write(reinterpret_cast<char *>(&u), sizeof(u));
          u_count_int++;
        }
    }
  if (is_two_boundaries)
    u_count_int +=  symbol_table.endo_nbr();
  for (j = 0; j < (int) symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  for (j = 0; j < (int) symbol_table.endo_nbr(); j++)
    SaveCode.write(reinterpret_cast<char *>(&j), sizeof(j));
  SaveCode.close();
}

void
ModelTree::writeLatexModelFile(const string &filename, ExprNodeOutputType output_type) const
{
  ofstream output;
  output.open(filename.c_str(), ios::out | ios::binary);
  if (!output.is_open())
    {
      cerr << "ERROR: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  output << "\\documentclass[10pt,a4paper]{article}" << endl
         << "\\usepackage[landscape]{geometry}" << endl
         << "\\usepackage{fullpage}" << endl
         << "\\begin{document}" << endl
         << "\\footnotesize" << endl;

  // Write model local variables
  for (map<int, expr_t>::const_iterator it = local_variables_table.begin();
       it != local_variables_table.end(); it++)
    {
      int id = it->first;
      expr_t value = it->second;

      output << "\\begin{equation*}" << endl
             << symbol_table.getName(id) << " = ";
      // Use an empty set for the temporary terms
      value->writeOutput(output, output_type);
      output << endl << "\\end{equation*}" << endl;
    }

  for (int eq = 0; eq < (int) equations.size(); eq++)
    {
      output << "\\begin{equation}" << endl
             << "% Equation " << eq+1 << endl;
      // Here it is necessary to cast to superclass ExprNode, otherwise the overloaded writeOutput() method is not found
      dynamic_cast<ExprNode *>(equations[eq])->writeOutput(output, output_type);
      output << endl << "\\end{equation}" << endl;
    }

  output << "\\end{document}" << endl;

  output.close();
}

void
ModelTree::addEquation(expr_t eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != NULL && beq->get_op_code() == oEqual);

  equations.push_back(beq);
}

void
ModelTree::addEquationTags(int i, const string &key, const string &value)
{
  equation_tags.push_back(make_pair(i, make_pair(key, value)));
}

void
ModelTree::addAuxEquation(expr_t eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);
  assert(beq != NULL && beq->get_op_code() == oEqual);

  aux_equations.push_back(beq);
}

void
ModelTree::addTrendVariables(vector<int> trend_vars, expr_t growth_factor) throw (TrendException)
{
  while (!trend_vars.empty())
    if (trend_symbols_map.find(trend_vars.back()) != trend_symbols_map.end())
      throw TrendException(symbol_table.getName(trend_vars.back()));
    else
      {
        trend_symbols_map[trend_vars.back()] = growth_factor;
        trend_vars.pop_back();
      }
}

void
ModelTree::addNonstationaryVariables(vector<int> nonstationary_vars, expr_t deflator) throw (TrendException)
{
  while (!nonstationary_vars.empty())
    if (nonstationary_symbols_map.find(nonstationary_vars.back()) != nonstationary_symbols_map.end())
      throw TrendException(symbol_table.getName(nonstationary_vars.back()));
    else
      {
        nonstationary_symbols_map[nonstationary_vars.back()] = deflator;
        nonstationary_vars.pop_back();
      }
}

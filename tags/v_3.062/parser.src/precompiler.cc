#include <map>
#include <stack>
#include <iostream>

#include "precompiler.hh"

extern "C" void start_skipping(void);
extern "C" void pop_flex_state(void);

static precompiler dyn_pc;

void precompiler::comp_define(std::string name, std::string integer)
{
  // go back to previous state before @define
  pop_flex_state();

  if (dyn_pc.macros.find(name) != dyn_pc.macros.end())
    {
      std::cout << "Warning: symbol " << name << " redefined" << "\n";
    }
  dyn_pc.macros[name]=std::atoi(integer.c_str());
}

void precompiler::comp_if(std::string name, std::string op, std::string integer)
{
  //go back to previous state before @if
  pop_flex_state();

  dyn_pc.precompiler_state.push(dyn_pc.current_state);
  dyn_pc.previous_skip_else.push(dyn_pc.skip_else);

  bool state;
  int value = atoi(integer.c_str());

  if ( op == "==" )
    {
      state = (dyn_pc.macros[name] == value);
    }
  else if ( op == "!=" )
    {
      state = (dyn_pc.macros[name] != value);
    }
  else if ( op == "<" )
    {
      state = (dyn_pc.macros[name] < value);
    }
  else if ( op == ">" )
    {
      state = (dyn_pc.macros[name] > value);
    }
  else if ( op == "<=" )
    {
      state = (dyn_pc.macros[name] <= value);
    }
  else if ( op == ">=" )
    {
      state = (dyn_pc.macros[name] >= value);
    }
  else
    {
      std::cout << "Precompiler: unkown operator";
    }

  dyn_pc.precompiler_state.push(dyn_pc.current_state);
  if ( dyn_pc.current_state == true && state == false)
    {
      start_skipping();
    }

  // it is possible to skip following else and elseif statements
  // if the previous state is already false (skipping is on)
  dyn_pc.skip_else = (dyn_pc.current_state) ? false : true;
  dyn_pc.current_state = state && dyn_pc.current_state;

}
  
void precompiler::comp_elseif(std::string name, std::string op, std::string integer)
{
  // go back to previous state before @elseif
  pop_flex_state();

  if (dyn_pc.skip_else)
    {
      return;
    }

  if (dyn_pc.current_state)
    {
      start_skipping();
      dyn_pc.current_state = false;
      dyn_pc.skip_else = true;
    }
  
  bool state;
  int value = atoi(integer.c_str());

  if ( op == "==" )
    {
      state = (dyn_pc.macros[name] == value);
    }
  else if ( op == "!=" )
    {
      state = (dyn_pc.macros[name] != value);
    }
  else if ( op == "<" )
    {
      state = (dyn_pc.macros[name] < value);
    }
  else if ( op == ">" )
    {
      state = (dyn_pc.macros[name] > value);
    }
  else if ( op == "<=" )
    {
      state = (dyn_pc.macros[name] <= value);
    }
  else if ( op == ">=" )
    {
      state = (dyn_pc.macros[name] >= value);
    }
  else
    {
      std::cout << "Precompiler: unkown operator";
    }

  if (state)
    {
      pop_flex_state();
      dyn_pc.current_state = true;
    }
}

void precompiler::comp_else(void)
{
  if (dyn_pc.skip_else)
    {
      return;
    }

  if (dyn_pc.current_state)
    {
      start_skipping();
      dyn_pc.current_state = false;
      dyn_pc.skip_else = true;
    }
  else  
    {
      pop_flex_state();
      dyn_pc.current_state = true;
    }
}  

void precompiler::comp_endif(void)
{
  if (!dyn_pc.current_state && dyn_pc.precompiler_state.top())
    {
      pop_flex_state();
    }

  dyn_pc.precompiler_state.pop();
  dyn_pc.current_state = dyn_pc.precompiler_state.top();

  dyn_pc.previous_skip_else.pop();
  dyn_pc.skip_else = dyn_pc.previous_skip_else.top();
}


extern "C" void compile_define(char* name, char* integer)
                                   { dyn_pc.comp_define(name, integer); }

extern "C" void compile_if(char* name, char* op, char* integer)
                                   { dyn_pc.comp_if(name,op,integer);}

extern "C" void compile_elseif(char* name, char* op, char* integer)
                                   { dyn_pc.comp_elseif(name,op,integer);}

extern "C" void compile_else(void) { dyn_pc.comp_else();}

extern "C" void compile_endif(void) { dyn_pc.comp_endif();}


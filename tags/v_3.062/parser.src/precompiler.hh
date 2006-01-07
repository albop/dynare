#include<string>

class precompiler {
  std::map<std::string,int> macros;
  std::stack<bool> precompiler_state;
  std::stack<bool> previous_skip_else;
  bool current_state;
  bool skip_else;
public:
  precompiler() 
  {
    current_state = true;
    skip_else = false;
    precompiler_state.push(current_state);
    previous_skip_else.push(skip_else);
  }

  void comp_define(std::string name, std::string integer);
  void comp_if(std::string name, std::string op, std::string integer);
  void comp_elseif(std::string name, std::string op, std::string integer);
  void comp_else(void);
  void comp_endif(void);
};

extern "C" void compile_define(char* name, char* integer);
extern "C" void compile_if(char* name, char* op, char* integer);
extern "C" void compile_elseif(char* name, char* op, char* integer);
extern "C" void compile_else(void);
extern "C" void compile_endif(void);



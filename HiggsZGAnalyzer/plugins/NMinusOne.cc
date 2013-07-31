#include "NMinusOne.h"


NMinusOne::NMinusOne(vector<void*> cuts, vector<void*> values, vector<string> ops){
  _cuts=cuts;
  _inputValues=values;
  _ops=ops;

  bool cutPass = false;

  //begin looping through cuts (only need to make sure this loop is N elements long, doesn't matter what you're looping over)
  for (unsigned int j = 0; j<cuts.size(); j++){
    
    //begin looping all N-1 cuts for the above cut
    for (unsigned int i =1; i<_cuts.size(); i++){
      cutPass = TestCut(_cuts[i],_inputValues[i],_ops[i]);
      if (!cutPass) _outputValues.push_back("fail!"); break;
    }
    if (cutPass) _outputValues.push_back(_inputValues[0]);

  //cycle the first element in all these vectors to the end to get ready for the next cut
  _cuts.push_back(_cuts[0]);
  _cuts.erase(_cuts.begin());
  _inputValues.push_back(_inputValues[0]);
  _inputValues.erase(_inputValues.begin());
  _ops.push_back(_ops[0]);
  _ops.erase(_ops.begin());

  }

}

bool NMinusOne::TestCut(void* cut, void* value, string op){
  bool _cutPass = false;

  if        (op == "less"){
    if (*value < *cut) _cutPass = true;
  }else if  (op == "greater"){
    if (*cut < *value) _cutPass = true;
  }else if  (op == "equal"){
    if (*value == *cut) _cutPass = true;
  }else if  (op == "!equal"){
    if (!(*value == *cut)) _cutPass = true;
  }else{
    cout<<"not a vaild operator: "<<op<<endl;
  }

  return _cutPass;
}
/*
vector<multiType> NMinusOne::GetOutput(){
  return _outputValues;
}
*/


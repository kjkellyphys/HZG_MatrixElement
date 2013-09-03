#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

int main()
{
  double DiscArray[200001];
  int TypeArray[200001];
  ifstream InputFile;
  InputFile.open ("Discriminant.txt");

  for (int i = 0; i < 200000; i++)
    {
      InputFile>>DiscArray[i];
      InputFile>>TypeArray[i];
    }

  InputFile.close();

  ofstream CutOutput;
  CutOutput.open ("DiscCuts.txt");

  double MaxSB, MaxSB_cut;
  int N_s, N_b;
  double Acc_s, Acc_b;
  double S, B;

  MaxSB = 0.; MaxSB_cut = 0.;

  for (double DiscCut = 0.; DiscCut < 1.0; DiscCut += 0.001){
    N_s = 0;
    N_b = 0;
    S = 0; B = 0;

    for (int i = 0; i < 200000; i++){
      if ((DiscArray[i] > DiscCut) && (TypeArray[i] == 0))
	N_s++;
      if ((DiscArray[i] > DiscCut) && (TypeArray[i] == 1))
	N_b++;
    }

    Acc_s = (double)N_s / 100000;
    Acc_b = (double)N_b / 100000;

    S = 0.0007091 * Acc_s;
    B = 22.89 * Acc_b;

    CutOutput << DiscCut << " " << N_s << " " << N_b << " " << (double)S/B << " " << sqrt(1/((double)S+B)) << " " << sqrt((double)(S+B)/S) << endl;

    if ((double) S/B > MaxSB){
      MaxSB = (double) S/B;
      MaxSB_cut = DiscCut;
    }
  }

  CutOutput << "Maximum S/B: " << MaxSB << " at Discriminant cut of " << MaxSB_cut << endl;
  CutOutput.close();
  return 0;
}

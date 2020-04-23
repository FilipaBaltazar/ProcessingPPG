import 'package:testing/testing.dart';

void main(List<String> arguments) {
  List<double> X = [0,1,2,3,4,5];
  print(X);
  List<double> Y = PPG.movingAverage(X,2);
  print(Y);
  X.add(6);
  print(X);
  PPG.addToMovingAverage(Y, X, 2);
  print(Y);
}

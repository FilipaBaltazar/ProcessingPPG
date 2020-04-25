import 'package:testing/testing.dart';
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';


void main(List<String> arguments) {
  var X = <double>[0,1,2,3,4,5];
  print(X);
  List<double> Y = [];
  PPG.completeMovingAverage(Y,X,3);
  print(Y);
  X.add(6);
  X.add(7);
  print(X);
  PPG.completeMovingAverage(Y, X, 3);
  print(Y);

  var N = 3.0;
  var fs = 128;
  var n = linspace(0, N, num: (N * fs).toInt(), endpoint: false);
  print(n);
  var f1 = 1; // 1Hz
  var sg1 = arraySin(arrayMultiplyToScalar(n, 2 * pi * f1));
  var f2 = 10; // 10Hz
  var sg2 = arraySin(arrayMultiplyToScalar(n, 2 * pi * f2));
  var sg = sg1 + sg2;
  // design a FIR filter low pass with cut frequency at 1Hz, the objective is
  // remove the 10Hz sine wave signal
  var nyq = 0.5 * fs; // nyquist frequency
  var fc = Array([0.5, 6]); // cut frequency 1Hz
  var normal_fc = fc / Array([nyq, nyq]); // frequency normalization for digital filters
  var numtaps = 30; // attenuation of the filter after cut frequency

  // generation of the filter coefficients
  var b = firwin(numtaps, normal_fc, pass_zero: false );
  print('FIR');
  print(b);

  //-------- Digital filter application -----------//
print('digital filter application');
// Apply the filter on the signal using lfilter function
// lfilter uses direct form II transposed, for FIR filter
// the a coefficient is 1.0
var sgFiltered = lfilter(b, Array([1.0]), sg);

print(sg); // show the original signal
print(sgFiltered); // show the filtered signal

//-------- peaks -----------//
print('peaks');
var pk = findPeaks(sgFiltered);
print(pk[0]);// print the indexes of the peaks found in the signal
print(pk[1]);// print the values of the peaks found in the signal



}

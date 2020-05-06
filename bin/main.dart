import 'package:testing/testing.dart';
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';


void main(List<String> arguments) {
  var X = <double>[0,1,2,3,4,5];
  print(X);
  
  var sampleInterval = Duration(milliseconds: 1000);
  var last = DateTime.now();

  var Y = <DateTime>[];

  PPG ppg = PPG(20);

  for (var i = 0; i < X.length; i++) {
    Y.add(last);
    ppg.add(last, X.elementAt(i));
    last = last.add(sampleInterval);
  }

  print(ppg.pulseRate);
}

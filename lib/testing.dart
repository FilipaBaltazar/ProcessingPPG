import 'package:smart_arrays_base/smart_arrays_base.dart';
import 'package:smart_arrays_numerics/smart_arrays_numerics.dart';
import 'package:smart_signal_processing/smart_signal_processing.dart';
import "dart:typed_data";
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';


class SensorValue {
  final DateTime time;
  final double value;

  SensorValue(this.time, this.value);
}

class PPG {
  List<double> valuesRaw = [];
  List<double> _valuesDC = [];
  List<double> _valuesAC = [];
  List<SensorValue> _sensorValuesAC = [];
  List<DateTime> times = [];
  int samplingRate;

  PPG(this.samplingRate);

  void add(time, value) {
    valuesRaw.add(value);
    times.add(time);
  }

  double get samplingInterval {
    return 1000 / samplingRate;
  }

  int get windowSize {
    return (60).round();
  }

  int get length {
    return valuesRaw.length;
  }

  List<double> get valuesDC {
    completeMovingAverage(_valuesDC, valuesRaw, windowSize);
    return _valuesDC;
  }

  List<double> get valuesAC {
    /* Subtraction of arrays */
    completeAddVV(_valuesAC, valuesRaw.sublist((windowSize/2).floor(),length-(windowSize/2).ceil()+1),valuesDC, type:'-');
    return _valuesAC;
  }

  List<DateTime> get timesAC {
    return times.sublist((windowSize/2).floor(),length-(windowSize/2).ceil()+1);
  }

  List<SensorValue> get sensorValuesAC {
    int i = _sensorValuesAC.length;
    while (i < valuesAC.length) {
      _sensorValuesAC.add(SensorValue(timesAC.elementAt(i), valuesAC.elementAt(i)));
      i++;
    }
    return _sensorValuesAC;
  }

  static void addToMovingAverage(List<double> Y, List<double> X, int N) {
    assert(Y.length == X.length - N);
    Y.add( Y.last - X.elementAt(X.length-N-1)/N + X.last/N );
  }

  static void completeMovingAverage(List<double> Y, List<double> X, int N) {
    var i = Y.length;
    if (i == 0) {
      var sumX = X.sublist(0, N).reduce((value, element) => value + element);
      Y.add(sumX/N);
      i++;
    }
    while (i < X.length-N+1) {
      Y.add( Y.elementAt(i-1) - X.elementAt(i-1)/N + X.elementAt(i+N-1)/N );
      i++;
    }
  }

  static List<double> movingAverage(List<double> X, int N){
    List<double> Y = [];
    var sumX = X.sublist(0, N).reduce((value, element) => value + element);
    Y.add(sumX/N);
    int i = 1;
    while (i < X.length-N+1) {
      Y.add( Y.elementAt(i-1) - X.elementAt(i-1)/N + X.elementAt(i+N-1)/N );
      i++;
    }
    return Y;
  }

  static void completeAddVV(List<double> Z, List<double> Y, List<double> X, {String type = '+'}){
    assert(Y.length == X.length);
    assert(type == '+' || type == '-');
    int i = Z.length;
    while (i < Y.length) {
      Z.add( type == '+' ? Y.elementAt(i)+X.elementAt(i) : Y.elementAt(i)-X.elementAt(i) );
    }
  }


 //static Array

}

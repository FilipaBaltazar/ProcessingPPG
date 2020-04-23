import 'package:smart_arrays_base/smart_arrays_base.dart';
import 'package:smart_arrays_numerics/smart_arrays_numerics.dart';
import 'package:smart_signal_processing/smart_signal_processing.dart';
import "dart:typed_data";

class SensorValue {
  final DateTime time;
  final double value;

  SensorValue(this.time, this.value);
}


class PPG {
  List<double> valuesRaw = [];
  List<DateTime> times = [];
  int samplingRate;

  PPG(this.samplingRate);

  void add(time, value) {
    this.valuesRaw.add(value);
    this.times.add(time);
  }

  double get samplingInterval {
    return 1000 / samplingRate;
  }

  int get windowSize {
    return (samplingRate*2).round();
  }

  int get length {
    return valuesRaw.length;
  }

  List<double> get valuesDC {
     return movingAverage(valuesRaw, windowSize);
  }

  List<double> get valuesAC {
    /* Subtraction of arrays */
    List<double> negDC = Numerics().dotCM(-1, [valuesDC]).elementAt(0);
    return Numerics().addVV(valuesRaw.sublist(windowSize), negDC);
  }

  List<DateTime> get timesAC {
    return times.sublist(windowSize);
  }

  List<SensorValue> get sensorValuesAC {
    List<SensorValue> result = [];
    int i = 0;
    while (i < valuesAC.length) {
      result.add(SensorValue(timesAC.elementAt(i), valuesAC.elementAt(i)));
      i++;
    }
    return result;
  }

  static void addToMovingAverage(List<double> Y, List<double> X, int N) {
    assert(Y.length == X.length - N);
    print(X.elementAt(X.length-N-1));
    Y.add( Y.last - X.elementAt(X.length-N-1)/N + X.last/N );
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

}
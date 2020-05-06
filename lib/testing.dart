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
  List<int> durations = [];
  DateTime start;
  int samplingRate;
  List<BasisFunction> _interpolation = [];
  List<double> _valuesInterp;
  List<double> _durationsInterp;
  List<double> _valuesFiltered;
  double _pulseRate;

  PPG(this.samplingRate);

  void add(DateTime time, double value) {
    valuesRaw.add(value);
    times.add(time);
    if (times.length == 1) {
      start = time;
    }
    durations.add(time.difference(start).inMilliseconds);
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

  List<double> get durationsInterp {
    if (_durationsInterp!=null) {
      return _durationsInterp;
    }
    else {
      _durationsInterp = List<double>.generate(
        (durations.last*samplingRate/1000).floor(), 
        (index) => index * 1/samplingRate*1000);
      return _durationsInterp;
    }
  }

  List<double> get valuesInterp {
    if (_valuesInterp!=null){
      return _valuesInterp;
    }
    else{
      _buildBasisFunctions();
      _valuesInterp = interpolate(durationsInterp);
      return _valuesInterp;
    }
  }

  List<double> get valuesFiltered {
    if (_valuesFiltered == null) {
      var nyq = 0.5 * samplingRate; // nyquist frequency
      var fc = Array([0.1, 3]); // cut frequency 1Hz
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
      _valuesFiltered = lfilter(b, Array([1.0]), Array(valuesInterp)).toList();
    }
    return _valuesFiltered;
  }

  double get pulseRate {
    if (_pulseRate == null) {
      print('peaks');
      var pk = findPeaks(valuesFiltered);
      var indices = pk[0];
      var timeSpan = durationsInterp.elementAt(indices.last) - durationsInterp.elementAt(indices.first);
      _pulseRate = (indices.length-1) / timeSpan / 1000 * 60;
    }
    return _pulseRate;
  }

  void _buildBasisFunctions() {
    _interpolation.add(
      BasisFunction(
        null, 
        durations.first, 
        durations.elementAt(1),
        valuesRaw.first
      )
    );
    for (var i=1; i<length-1; i++) {
      _interpolation.add(
        BasisFunction(
          durations.elementAt(i-1), 
          durations.elementAt(i), 
          durations.elementAt(i+1),
          valuesRaw.elementAt(i)
        )
      );
    }
    _interpolation.add(
      BasisFunction(
        durations.elementAt(length-2), 
        durations.last, 
        null,
        valuesRaw.last
      )
    );
  }

  List<num> interpolate(List<num> times) {
    assert(times.first >= _interpolation.first.midTime);
    assert(times.last <= _interpolation.last.finalTime);
    var interpIndex = 0;
    var timesIndex = 0;
    var result = [];
    while (timesIndex < times.length) {
      if (times.elementAt(timesIndex) <= _interpolation.elementAt(interpIndex).finalTime) {
        result.add(_interpolation.elementAt(interpIndex).midTimeValue+_interpolation.elementAt(interpIndex+1).midTimeValue);
        timesIndex++;
      }
      else {
        interpIndex++;
      }
    }
    return result;
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

class BasisFunction {
  num startTime;
  num midTime;
  num finalTime;
  num midTimeValue;

  BasisFunction(this.startTime, this.midTime, this.finalTime, this.midTimeValue);

  num at(num t) {
    if (t < startTime) {
      return 0;
    }
    else if (t < midTime) {
      return midTimeValue * (t-startTime) / (midTime-startTime);
    }
    else if (t < finalTime) {
      return midTimeValue * (finalTime-t) / (finalTime-midTime);
    }
    else {
      return 0;
    }
  }
}
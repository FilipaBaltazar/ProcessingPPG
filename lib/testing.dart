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

class Oximetry {
  PPG signalRed;
  PPG signalGreen;
  PPG signalBlue;
  // wavelengths in nm
  int lambdaRed = 640;
  int lambdaGreen = 520;
  int lambdaBlue = 450;
  // coefficients of extinction in L / mmol / cm
  double eHbOxyRed = 0.442;
  double eHbOxyGreen = 2.42024;
  double eHbOxyBlue = 6.2816;
  double eHbRed = 4.3542;
  double eHbGreen = 3.1589;
  double eHbBlue = 10.3292;

  Oximetry(this.signalRed, this.signalBlue, [this.signalGreen]);

  /// Returns the SpO2 value of [this].
  double get value {
    var paramsRed = signalParams(signalRed);
    var slopeRed = mean(paramsRed[0]);
    var peakRed = mean(paramsRed[1]);
    var paramsBlue = signalParams(signalBlue);
    var slopeBlue = mean(paramsBlue[0]);
    var peakBlue = mean(paramsBlue[1]);
    return 100
     * ( eHbRed*sqrt(slopeBlue*peakBlue) - eHbBlue*sqrt(slopeRed*peakRed) )
     / ( sqrt(slopeBlue*peakBlue)*(eHbRed - eHbOxyRed) - sqrt(slopeRed*peakRed)*(eHbBlue - eHbOxyBlue) );
  }

  /// Returns the list of slopes and peak to peak values of [signal], in that order.
  static List<Array> signalParams(PPG signal) {
    var peak2peaks = Array.empty();
    var slopes = Array.empty();
    for (var i=1; i<signal.peaks[0].length; i++) {
      // index of the prev positive peak
      var indHighPrev = signal.peaks[0][i-1];
      // index of positive peak
      var indHigh = signal.peaks[0][i];
      // index of negative peak
      var indLow = signal.negativePeaks[0].indexWhere((var element) => indHighPrev < element && element < indHigh);
      if (indLow != -1) {
        // timepoint of positive peak
        var timeHigh = signal.durationsInterp[indHigh];
        // timepoint of negative peak
        var timeLow = signal.durationsInterp[indLow];
        // value of positive peak
        var valueHigh = signal.valuesLog[indHigh];
        // value of negative peak
        var valueLow = signal.valuesLog[indLow];
        var peak2peak = valueHigh - valueLow;
        var slope = peak2peak / (timeHigh - timeLow);
        peak2peaks.add(peak2peak);
        slopes.add(slope);
      }
    }
    return [slopes,peak2peaks];
  }
}

class PPG {
  Array valuesRaw = Array.empty();
  Array _valuesDC = Array.empty();
  Array _valuesAC = Array.empty();
  List<SensorValue> _sensorValuesAC = [];
  List<DateTime> times = [];
  Array durations = Array.empty();
  DateTime start;
  int samplingRate;
  List<BasisFunction> _interpolation = [];
  Array _valuesInterp;
  Array _durationsInterp;
  Array _valuesLog;
  double _pulseRate;
  List<dynamic> _peaks;
  List<dynamic> _negativePeaks;

  PPG(this.samplingRate);

  void add(DateTime time, double value) {
    valuesRaw.add(value);
    times.add(time);
    if (times.length == 1) {
      start = time;
    }
    durations.add(time.difference(start).inMilliseconds.toDouble());
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

  Array get valuesDC {
    completeMovingAverage(_valuesDC, valuesRaw, windowSize);
    return _valuesDC;
  }

  Array get valuesAC {
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

  List<BasisFunction> get interpolation {
    return _interpolation;
  }

  Array get durationsInterp {
    if (_durationsInterp!=null) {
      return _durationsInterp;
    }
    else {
      _durationsInterp = Array(List<double>.generate(
        (durations.last*samplingRate/1000).floor()+1, 
        (index) => index * 1/samplingRate*1000));
      return _durationsInterp;
    }
  }

  Array get valuesInterp {
    if (_valuesInterp!=null){
      return _valuesInterp;
    }
    else{
      _buildBasisFunctions();
      _valuesInterp = interpolate(durationsInterp);
      return _valuesInterp;
    }
  }

  Array get valuesLog {
    if (_valuesLog == null) {
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
      _valuesLog = arrayLog(lfilter(b, Array([1.0]), valuesInterp));
    }
    return _valuesLog;
  }

  List<dynamic> get peaks {
    _peaks ??= findPeaks(valuesLog);
    return _peaks;
  }

  List<dynamic> get negativePeaks {
    _negativePeaks ??= findPeaks(Array.fixed(valuesLog.length)-valuesLog);
    return _negativePeaks;
  }

  double get pulseRate {
    if (_pulseRate == null) {
      print('peaks');
      var indices = peaks[0];
      print('Filtered values');
      print(valuesLog);
      // print(durationsInterp.last);
      // print(durations.last);
      var timeSpan = durationsInterp[(indices.last.round())] - durationsInterp[(indices.first.round())];
      _pulseRate = ((indices.length-1) / (timeSpan / 1000 ))* 60;
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

  Array interpolate(Array times) {
    assert(times.first >= _interpolation.first.midTime);
    assert(times.last <= _interpolation.last.midTime);
    var interpIndex = 0;
    var timesIndex = 0;
    var result = Array.empty();
    while (timesIndex < times.length) {
      if (times[timesIndex] <= _interpolation.elementAt(interpIndex).finalTime) {
        result.add(
          _interpolation.elementAt(interpIndex).at(times[timesIndex])
          +_interpolation.elementAt(interpIndex+1).at(times[timesIndex])
        );
        timesIndex++;
      }
      else {
        interpIndex++;
      }
    }
    return result;
  }

  static void addToMovingAverage(Array Y, Array X, int N) {
    assert(Y.length == X.length - N);
    Y.add( Y.last - X.elementAt(X.length-N-1)/N + X.last/N );
  }

  static void completeMovingAverage(Array Y, Array X, int N) {
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

  static Array movingAverage(Array X, int N){
    var Y = <double>[];
    var sumX = X.sublist(0, N).reduce((value, element) => value + element);
    Y.add(sumX/N);
    var i = 1;
    while (i < X.length-N+1) {
      Y.add( Y.elementAt(i-1) - X.elementAt(i-1)/N + X.elementAt(i+N-1)/N );
      i++;
    }
    return Y;
  }

  static void completeAddVV(Array Z, Array Y, Array X, {String type = '+'}){
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
    if (startTime != null && t < startTime) {
      return 0.0;
    }
    else if (t < midTime) {
      return midTimeValue * (t-startTime) / (midTime-startTime);
    }
    else if (t == midTime) {
      return midTimeValue;
    }
    else if (finalTime != null && t < finalTime) {
      return midTimeValue * (finalTime-t) / (finalTime-midTime);
    }
    else {
      return 0.0;
    }
  }

  @override
  String toString() {
    return 'Start time:'+startTime.toString()+' Mid time: '+midTime.toString()+' finalTime: '+finalTime.toString()+' Value: '+midTimeValue.toString();
  }

}
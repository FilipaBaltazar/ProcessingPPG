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
  /// Wavelength assumed for the red channel, in nm.
  static int lambdaRed = 640;

  /// Wavelength assumed for the green channel, in nm.
  static int lambdaGreen = 520;

  /// Wavelength assumed for the blue channel, in nm.
  static int lambdaBlue = 450;

  /// Extinction coefficient of HbO at the red wavelength, in L / mmol / cm.
  static double eHbOxyRed = 0.442;

  /// Extinction coefficient of HbO at the green wavelength, in L / mmol / cm.
  static double eHbOxyGreen = 2.42024;

  /// Extinction coefficient of HbO at the blue wavelength, in L / mmol / cm.
  static double eHbOxyBlue = 6.2816;

  /// Extinction coefficient of Hb at the red wavelength, in L / mmol / cm.
  static double eHbRed = 4.3542;

  /// Extinction coefficient of Hb at the green wavelength, in L / mmol / cm.
  static double eHbGreen = 3.1589;

  /// Extinction coefficient of Hb at the blue wavelength, in L / mmol / cm.
  static double eHbBlue = 10.3292;

  /// Returns the SpO2 value of calculated from the [PPG] signals `signalRed` and `signalBlue`.
  ///
  /// The value is calculated using the method described by Reddy et al. (2009).
  static double value(PPG signalRed, PPG signalBlue) {
    var paramsRed = signalParams(signalRed);
    var slopeRed = mean(paramsRed[0]);
    var peakRed = mean(paramsRed[1]);
    var paramsBlue = signalParams(signalBlue);
    var slopeBlue = mean(paramsBlue[0]);
    var peakBlue = mean(paramsBlue[1]);
    return 100 *
        (eHbRed * sqrt(slopeBlue * peakBlue) -
            eHbBlue * sqrt(slopeRed * peakRed)) /
        (sqrt(slopeBlue * peakBlue) * (eHbRed - eHbOxyRed) -
            sqrt(slopeRed * peakRed) * (eHbBlue - eHbOxyBlue));
  }

  /// Returns the list of slopes and peak-to-peak values of `signal`, in that order.
  /// 
  /// For each pair of consecutive positive peaks, 
  /// the function checks whether there is a negative peak between them.
  /// If there is, then the slope and peak-to-peak value are calculated
  /// and added to their respective lists.
  /// In the case that there is more than one negative peak, the one with lowest index is considered.
  /// The function terminates when all positive peaks have been checked.
  static List<Array> signalParams(PPG signal) {
    var peak2peaks = Array.empty();
    var slopes = Array.empty();
    for (var i = 1; i < signal.peaks[0].length; i++) {
      // index of the prev positive peak
      var indHighPrev = signal.peaks[0][i - 1];
      // index of positive peak
      var indHigh = signal.peaks[0][i];
      // index of negative peak
      var indLow = signal.negativePeaks[0].indexWhere(
          (var element) => indHighPrev < element && element < indHigh);
      if (indLow != -1) {
        // timepoint of positive peak
        var timeHigh = signal.millisInterp[indHigh];
        // timepoint of negative peak
        var timeLow = signal.millisInterp[indLow];
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
    return [slopes, peak2peaks];
  }
}

class PPG {
  /// Raw values of each data point in the PPG signal.
  Array valuesRaw = Array.empty();

  Array _valuesDC = Array.empty();
  Array _valuesAC = Array.empty();
  List<SensorValue> _sensorValuesAC = [];

  /// [DateTime] values of each data point in the PPG signal.
  List<DateTime> dates = [];

  /// First value in [dates].
  DateTime start;

  /// Timespan since [start] of each data point in the PPG signal, in milliseconds.
  Array millis = Array.empty();

  /// Sampling rate used to interpolate the PPG signal.
  int samplingRate;

  final List<BasisFunction> _interpolation = [];
  Array _valuesInterp;
  Array _millisInterp;
  Array _valuesLog;
  double _pulseRate;
  List<dynamic> _peaks;
  List<dynamic> _negativePeaks;

  /// High frequency cut-off of the band-pass filter applied to
  /// [valuesInterp] before obtaining [valuesLog].
  double frequencyHigh = 3;

  /// Low frequency cut-off of the band-pass filter applied to
  /// [valuesInterp] before obtaining [valuesLog].
  double frequencyLow = 0.1;

  PPG(this.samplingRate);

  void add(DateTime time, double value) {
    valuesRaw.add(value);
    dates.add(time);
    if (dates.length == 1) {
      start = time;
    }
    millis.add(time.difference(start).inMilliseconds.toDouble());
  }

  double get samplingInterval {
    return 1000 / samplingRate;
  }

  int get windowSize {
    return (60).round();
  }

  /// Length of [valuesRaw].
  int get length {
    return valuesRaw.length;
  }

  Array get valuesDC {
    completeMovingAverage(_valuesDC, valuesRaw, windowSize);
    return _valuesDC;
  }

  Array get valuesAC {
    /* Subtraction of arrays */
    completeAddVV(
        _valuesAC,
        valuesRaw.sublist(
            (windowSize / 2).floor(), length - (windowSize / 2).ceil() + 1),
        valuesDC,
        type: '-');
    return _valuesAC;
  }

  List<DateTime> get timesAC {
    return dates.sublist(
        (windowSize / 2).floor(), length - (windowSize / 2).ceil() + 1);
  }

  List<SensorValue> get sensorValuesAC {
    var i = _sensorValuesAC.length;
    while (i < valuesAC.length) {
      _sensorValuesAC
          .add(SensorValue(timesAC.elementAt(i), valuesAC.elementAt(i)));
      i++;
    }
    return _sensorValuesAC;
  }

  /// Returns the list of linear basis functions used to interpolate [valuesRaw].
  List<BasisFunction> get interpolation {
    if (valuesRaw.length > 2 && _interpolation.length < valuesRaw.length) {
      _completeBasisFunctions();
    }
    return _interpolation;
  }

  /// Returns the timespan since [start] of each data point in [valuesInterp], in milliseconds.
  Array get millisInterp {
    if (_millisInterp == null) {
      _millisInterp = Array(List<double>.generate(
          (millis.last * samplingRate / 1000).floor() + 1,
          (index) => index * 1 / samplingRate * 1000));
    }
    return _millisInterp;
  }

  /// Returns the interpolated values obtained from [valuesRaw] at a sampling rate of [samplingRate].
  Array get valuesInterp {
    if (_valuesInterp == null) {
      _valuesInterp = interpolate(millisInterp);
    }
    return _valuesInterp;
  }

  /// Returns the logarithm of [valuesInterp], applied after filtering.
  Array get valuesLog {
    if (_valuesLog == null) {
      var frequencyNyquist = 0.5 * samplingRate;
      var bandPassWindow = Array([frequencyLow, frequencyHigh]);
      var normalBandPassWindow =
          bandPassWindow / Array([frequencyNyquist, frequencyNyquist]);
      var filterOrder = 30;
      var filterCoeffs =
          firwin(filterOrder, normalBandPassWindow, pass_zero: false);
      _valuesLog = arrayLog(lfilter(filterCoeffs, Array([1.0]), valuesInterp));
    }
    return _valuesLog;
  }

  /// Returns the indices and the values of the positive peaks of [valuesLog],
  /// in that order.
  /// 
  /// A value is considered a positive peak if both 
  /// the previous and next value are smaller or equal to it.
  List<dynamic> get peaks {
    if (_peaks == null) {
      _peaks = findPeaks(valuesLog);
      _peaks[0] = _peaks[0].map((i) => i.round()).toList();
    }
    return _peaks;
    
  }

  /// Returns the indices and the values of the negative peaks of [valuesLog],
  /// in that order.
  /// 
  /// A value is considered a negative peak if both 
  /// the previous and next value are greater or equal to it.
  List<dynamic> get negativePeaks {
    if (_negativePeaks == null) {
      var zeroes = Array(List<double>.filled(valuesLog.length, 0, growable: true));
      _negativePeaks = findPeaks(zeroes - valuesLog);
      _negativePeaks[0] = _negativePeaks[0].map((i) => i.round()).toList();
    }
    return _negativePeaks;
  }

  /// Returns the pulse rate calculated from [valuesLog].
  double get pulseRate {
    if (_pulseRate == null) {
      // print('peaks');
      var indices = peaks[0];
      // print('Filtered values');
      // print(valuesLog);
      // print(durationsInterp.last);
      // print(durations.last);
      var timeSpan = millisInterp[(indices.last.round())] -
          millisInterp[(indices.first.round())];
      _pulseRate = ((indices.length - 1) / (timeSpan / 1000)) * 60;
    }
    return _pulseRate;
  }

  /// Fills in [_interpolation] with [BasisFunction] objects until its length equals that of [valuesRaw].
  /// 
  /// This is only performed if the [valuesRaw] has more than two elements, since a basis function
  /// on its own is useless.
  void _completeBasisFunctions() {
    if (valuesRaw.length < 3) {
      return;
    }
    if (_interpolation.isEmpty) {
      _interpolation.add(BasisFunction(
          null, millis[0], millis[1], valuesRaw[0]));
      _interpolation.add(BasisFunction(
        millis[0], millis[1], null, valuesRaw[1]));
    }
    if (_interpolation.length < length) {
      _interpolation.removeLast();
      for (var i = interpolation.length; i < length - 1; i++) {
        _interpolation.add(BasisFunction(
            millis[i - 1],
            millis[i],
            millis[i + 1],
            valuesRaw[i]));
      }
      _interpolation.add(BasisFunction(
          millis[length - 2], millis.last, null, valuesRaw.last)); 
    }
  }

  Array interpolate(Array times) {
    var result = Array.empty();
    if (interpolation.isEmpty) return result;
    assert(times.first >= interpolation.first.midTime);
    assert(times.last <= interpolation.last.midTime);
    var indexInterp = 0;
    var indexTime = 0;
    while (indexTime < times.length) {
      if (times[indexTime] <= interpolation[indexInterp].finalTime) {
        result.add(interpolation[indexInterp].at(times[indexTime]) +
            interpolation[indexInterp + 1].at(times[indexTime]));
        indexTime++;
      } else {
        indexInterp++;
      }
    }
    return result;
  }

  static void addToMovingAverage(Array Y, Array X, int N) {
    assert(Y.length == X.length - N);
    Y.add(Y.last - X.elementAt(X.length - N - 1) / N + X.last / N);
  }

  static void completeMovingAverage(Array Y, Array X, int N) {
    var i = Y.length;
    if (i == 0) {
      var sumX = X.sublist(0, N).reduce((value, element) => value + element);
      Y.add(sumX / N);
      i++;
    }
    while (i < X.length - N + 1) {
      Y.add(Y.elementAt(i - 1) -
          X.elementAt(i - 1) / N +
          X.elementAt(i + N - 1) / N);
      i++;
    }
  }

  static Array movingAverage(Array X, int N) {
    var Y = <double>[];
    var sumX = X.sublist(0, N).reduce((value, element) => value + element);
    Y.add(sumX / N);
    var i = 1;
    while (i < X.length - N + 1) {
      Y.add(Y.elementAt(i - 1) -
          X.elementAt(i - 1) / N +
          X.elementAt(i + N - 1) / N);
      i++;
    }
    return Y;
  }

  static void completeAddVV(Array Z, Array Y, Array X, {String type = '+'}) {
    assert(Y.length == X.length);
    assert(type == '+' || type == '-');
    int i = Z.length;
    while (i < Y.length) {
      Z.add(type == '+'
          ? Y.elementAt(i) + X.elementAt(i)
          : Y.elementAt(i) - X.elementAt(i));
    }
  }

  //static Array

}

class BasisFunction {
  num startTime;
  num midTime;
  num finalTime;
  num midTimeValue;

  BasisFunction(
      this.startTime, this.midTime, this.finalTime, this.midTimeValue);

  /// Returns the value of the interpolation at time point `t`.
  num at(num t) {
    if (startTime != null && t < startTime) {
      return 0.0;
    } else if (t < midTime) {
      return midTimeValue * (t - startTime) / (midTime - startTime);
    } else if (t == midTime) {
      return midTimeValue;
    } else if (finalTime != null && t < finalTime) {
      return midTimeValue * (finalTime - t) / (finalTime - midTime);
    } else {
      return 0.0;
    }
  }

  /// Returns the value of the interpolation at time point `t`.
  num operator [](num t) => at(t);

  @override
  String toString() {
    return 'Start time:' +
        startTime.toString() +
        ' Mid time: ' +
        midTime.toString() +
        ' finalTime: ' +
        finalTime.toString() +
        ' Value: ' +
        midTimeValue.toString();
  }
}

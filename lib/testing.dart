import 'package:smart_arrays_base/smart_arrays_base.dart';
import 'package:smart_arrays_numerics/smart_arrays_numerics.dart';
import 'package:smart_signal_processing/smart_signal_processing.dart';
import "dart:typed_data";
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';
import 'dart:math';
import 'package:linalg/linalg.dart';

class SensorValue {
  final DateTime time;
  final double value;

  SensorValue(this.time, this.value);
}

class Oximetry {
  /// Wavelength assumed for the red channel, in nm.
  static int lambdaRed = 668; //Used to be 612

  /// Wavelength assumed for the green channel, in nm.
  static int lambdaGreen = 550;

  /// Wavelength assumed for the blue channel, in nm.
  static int lambdaBlue = 468;

  /// Extinction coefficient of HbO at the red wavelength, in L / mmol / cm.
  static double eHbOxyRed =298;
	; 	

  /// Extinction coefficient of HbO at the green wavelength, in L / mmol / cm.
  static double eHbOxyGreen = 43016;

  /// Extinction coefficient of HbO at the blue wavelength, in L / mmol / cm.
  static double eHbOxyBlue = 34870.8		;

  /// Extinction coefficient of Hb at the red wavelength, in L / mmol / cm.
  static double eHbRed = 2881.4	;

  /// Extinction coefficient of Hb at the green wavelength, in L / mmol / cm.
  static double eHbGreen = 53412;

  /// Extinction coefficient of Hb at the blue wavelength, in L / mmol / cm.
  static double eHbBlue =	17025.6	;

  /// Returns the SpO2 value of calculated from the [PPG] signals `signalRed` and `signalBlue`.
  ///
  /// The value is calculated using the method described by Reddy et al. (2009).
  static double value(PPG signalRed, PPG signalBlue) {
    var paramsRed = signalParams(signalRed, blueQ: false);
    var slopeRed = median(paramsRed[0]);
    var peakRed = median(paramsRed[1]);
    var paramsBlue = signalParams(signalBlue);
    var slopeBlue = median(paramsBlue[0]);
    var peakBlue = median(paramsBlue[1]);
    var ratio = (eHbRed * sqrt(slopeBlue * peakBlue) - eHbBlue * sqrt(slopeRed * peakRed)) /
        (eHbOxyBlue * sqrt(slopeRed * peakRed) - eHbOxyRed * sqrt(slopeBlue * peakBlue));
    // var ratio = (eHbRed * (peakBlue) - eHbBlue * (peakRed)) /
    //     (eHbOxyBlue * (peakRed) - eHbOxyRed * (peakBlue));
    return (100 * (ratio / (ratio + 1 )));
  }

  /// Returns the list of slopes and peak-to-peak values of `signal`, in that order.
  /// 
  /// For each pair of consecutive positive peaks, 
  /// the function checks whether there is a negative peak between them.
  /// If there is, then the slope and peak-to-peak value are calculated
  /// and added to their respective lists.
  /// In the case that there is more than one negative peak, the one with lowest index is considered.
  /// The function terminates when all positive peaks have been checked.
  static List<Array> signalParams(PPG signal, {bool blueQ = true}) {
    var peak2peaks = Array.empty();
    var slopes = Array.empty();
    for (var i = 1; i < signal.peaks[0].length; i++) {
      // index of the prev positive peak
      var indHighPrev = signal.peaks[0][i - 1];
      // index of positive peak
      var indHigh = signal.peaks[0][i];
      // index of negative peak
      var indLow = signal.valleys[0].firstWhere(
          (var element) => indHighPrev < element && element < indHigh,
          orElse: () => -1);
      if (indLow != -1) {
        // timepoint of previous positive peak
        var timeHighPrev = signal.millisInterp[indHighPrev];
        // timepoint of positive peak
        var timeHigh = signal.millisInterp[indHigh];
        // timepoint of negative peak
        var timeLow = signal.millisInterp[indLow];
        // value of previous positive peak
        var valueHighPrev = signal.valuesLog[indHighPrev];
        // value of positive peak
        var valueHigh = signal.valuesLog[indHigh];
        // value of negative peak
        var valueLow = signal.valuesLog[indLow];
        // var peak2peak = blueQ ? valueHigh - valueLow : valueHighPrev - valueLow;
        // var deltaT = blueQ ? timeHigh - timeLow : timeLow - timeHighPrev;
        var peak2peak =  valueHighPrev - valueLow;
        var deltaT =  timeLow - timeHighPrev;
        var slope = peak2peak / deltaT * 1000;
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
  List<SensorValue> _sensorValuesRaw = [];
  List<SensorValue> _sensorValuesLog = [];
  List<SensorValue> _sensorValuesInterp = [];

  /// [DateTime] values of each data point in the PPG signal.
  List<DateTime> dates = [];

  /// First value in [dates].
  DateTime start;

  /// Timespan since [start] of each data point in the PPG signal, in milliseconds.
  Array millis = Array.empty();

  /// Sampling rate used to interpolate the PPG signal.
  int samplingRate;

  final List<BasisFunction> _interpolation = [];
  final List<BasisFunction> _interpolationPeaks = [];
  final List<BasisFunction> _interpolationValleys = [];
  Array _valuesInterp = Array.empty();
  Array _millisInterp = Array.empty();
  Array _valuesFiltered = Array.empty();
  double _pulseRate;
  List<dynamic> _peaks = [];
  List<dynamic> _valleys = [];


  /// High frequency cut-off of the band-pass filter applied to
  /// [valuesInterp] before obtaining [valuesFiltered].
  double frequencyHigh = 2;

  /// Low frequency cut-off of the band-pass filter applied to
  /// [valuesInterp] before obtaining [valuesFiltered].
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

  

  List<SensorValue> get sensorValuesRaw {
    var i = _sensorValuesRaw.length;
    while (i < valuesRaw.length) {
      _sensorValuesRaw
          .add(SensorValue(dates[i], valuesRaw[i]));
      i++;
    }
    return _sensorValuesRaw;
  }

  List<SensorValue> get sensorValuesLog {
    var i = _sensorValuesLog.length;
    while (i < valuesFiltered.length) {
      var date = start.add(Duration(microseconds: (millisInterp[i]*1000).round()));
      _sensorValuesLog
          .add(SensorValue(date, valuesFiltered[i]));
      i++;
    }
    return _sensorValuesLog;
  }

  List<SensorValue> get sensorValuesInterp {
    var i = _sensorValuesInterp.length;
    while (i < valuesInterp.length) {
      var date = start.add(Duration(microseconds: (millisInterp[i]*1000).round()));
      _sensorValuesInterp
          .add(SensorValue(date, valuesInterp[i]));
      i++;
    }
    return _sensorValuesInterp;
  }

  /// Returns the list of linear basis functions used to interpolate [valuesRaw].
  List<BasisFunction> get interpolation {
    if (valuesRaw.length > 2 && _interpolation.length < valuesRaw.length) {
      _fillBasisFunctionsPPG();
    }
    return _interpolation;
  }

  List<BasisFunction> get interpolationPeaks {
    if (peaks.length > 2 && _interpolationPeaks.length < peaks.length) {
      _fillBasisFunctionsEnvelope();
    }
    return _interpolationPeaks;
  }

  List<BasisFunction> get interpolationValleys {
    if (valleys.length > 2 && _interpolationValleys.length < valleys.length) {
      _fillBasisFunctionsEnvelope();
    }
    return _interpolationValleys;
  }

  /// Returns the timespan since [start] of each data point in [valuesInterp], in milliseconds.
  Array get millisInterp {
    if (true) {
      _millisInterp = Array(List<double>.generate(
          (millis.last * samplingRate / 1000).floor() + 1,
          (index) => index * 1 / samplingRate * 1000));
    }
    return _millisInterp;
  }

  /// Returns the interpolated values obtained from [valuesRaw] at a sampling rate of [samplingRate].
  Array get valuesInterp {
    if (true) {
      _valuesInterp = interpolatePPG();
    }
    return _valuesInterp;
  }

  /// Returns a filtered version of [valuesInterp].
  Array get valuesFiltered {
    if (true) {
      var frequencyNyquist = 0.5 * samplingRate;
      var bandPassWindow = Array([frequencyLow, frequencyHigh]);
      var normalBandPassWindow =
          bandPassWindow / Array([frequencyNyquist, frequencyNyquist]);
      var filterOrder = 30;
      var filterCoeffs =
          firwin(filterOrder, normalBandPassWindow, pass_zero: false);
      _valuesFiltered = (lfilter(filterCoeffs, Array([1.0]), valuesInterp));
    }
    return _valuesFiltered;
  }

  /// Returns the natural logarithm of [valuesFiltered].
  Array get valuesLog {
    return arrayLog(valuesFiltered);
  }

  /// Returns the mean envelope of [valuesFiltered].
  Array get valuesMeanEnvelope {
    var _valuesUpper = interpolate(interpolationPeaks, millisInterp);
    var _valuesLower = interpolate(interpolationValleys, millisInterp);
    var result = (_valuesUpper + _valuesLower);
    for (var i=0; i<result.length; i++) {
      result[i] = result[i]/2;
    }
    return result;
  }

  /// Returns the difference between [valuesFiltered] and [valuesMeanEnvelope].
  Array get valuesProcessed {
    return valuesFiltered - valuesMeanEnvelope;
  }

  /// Returns the indices and the values of the positive peaks of [valuesFiltered],
  /// in that order.
  /// 
  /// A value is considered a positive peak if both 
  /// the previous and next value are smaller or equal to it.
  List<dynamic> get peaks {
    if (true) {
      _peaks = findPeaks(valuesFiltered);
      _peaks[0] = _peaks[0].map((i) => i.round()).toList();
      var millis = Array.empty();
      _peaks[0].forEach((i) => {millis.add(millisInterp[i])});
      _peaks.length < 3 ? _peaks.add(millis) : _peaks[3] = millis;
    }
    return _peaks;
    
  }

  /// Returns the indices and the values of the valleys of [valuesFiltered],
  /// in that order.
  /// 
  /// A value is considered a negative peak if both 
  /// the previous and next value are greater or equal to it.
  List<dynamic> get valleys {
    if (true) {
      _valleys = findValleys(valuesFiltered);
      _valleys[0] = _valleys[0].map((i) => i.round()).toList();
      var millis = Array.empty();
      _valleys[0].forEach((i) => {millis.add(millisInterp[i])});
      _valleys.length < 3 ? _valleys.add(millis) : _valleys[3] = millis;
    }
    return _valleys;
  }

  /// Returns the pulse rate calculated from [valuesFiltered].
  double get pulseRate {
    if (true) {
      // print('peaks');
      //var indices = peaks[0];
      // print('Filtered values');
      // print(valuesLog);
      // print(durationsInterp.last);
      // print(durations.last);
      var timeSpan = peaks[2].last-peaks[2].first;
      _pulseRate = ((peaks[2].length - 1) / (timeSpan / 1000)) * 60;
      print('Peak count pulse rate');
      print(_pulseRate);

      var intervalsRR = arrayDiff(peaks[2]);
      var ratesRR = Array.empty();
      intervalsRR.forEach((i) => {ratesRR.add(1/i*1000*60)});
      _pulseRate = mean(ratesRR);
      print('Mean RR inverse pulse rate');
      print(_pulseRate);

      // var fftinterp = arrayComplexAbs(rfft(valuesInterp,n:valuesInterp.length));
      // //var fftinterp = FFT.transform(Float64List.fromList(valuesInterp.toList()), Float64List.fromList([]));
      // print('fft interp');
      // print(fftinterp);
      // var fftfiltered = arrayComplexAbs(rfft(valuesFiltered,n:valuesFiltered.length));
      // print('fft filtered');
      // print(fftfiltered);
      // var fft = arrayComplexAbs(rfft(valuesProcessed,n:valuesProcessed.length));
      // print('fft processed');
      // print(fft);
      // var frequency = fftFreq(valuesInterp.length, d: 1/samplingRate/60);
      // var index = arrayArgMax(fft.getRangeArray(10,fft.length));
      // _pulseRate = frequency[index];
      // print('FFT pulse rate');
      // print(_pulseRate);

      // print('FFT frequency domain (BPM)');
      // print(frequency);
    

    }
    return _pulseRate;
  }

  /// Fills in [_interpolation] with [BasisFunction] objects until its length equals that of [valuesRaw].
  /// 
  /// This is only performed if the [valuesRaw] has more than two elements, since a basis function
  /// on its own is useless.
  void _fillBasisFunctionsPPG() {
    _fillBasisFunctions(_interpolation, millis, valuesRaw);
  }

  void _fillBasisFunctionsEnvelope(){
      var times = Array(<double>[0.0] + peaks[2].toList() + <double>[millisInterp.last]);
      var values = Array(<double>[valuesFiltered.first] + peaks[1].toList() + <double>[valuesFiltered.last]);
      _fillBasisFunctions(_interpolationPeaks, times, values);
      times = Array(<double>[0.0] + valleys[2].toList() + <double>[millisInterp.last]);
      values = Array(<double>[valuesFiltered.first] + valleys[1].toList() + <double>[valuesFiltered.last]);
      _fillBasisFunctions(_interpolationValleys, times, values);
  }

  static void _fillBasisFunctions(List<BasisFunction> basisFunctions, Array times, Array values) {
    if (values.length < 3) {
      return;
    }
    if (basisFunctions.isEmpty) {
      basisFunctions.add(BasisFunction(
          null, times[0], times[1], values[0]));
      basisFunctions.add(BasisFunction(
        times[0], times[1], null, values[1]));
    }
    if (basisFunctions.length < values.length) {
      basisFunctions.removeLast();
      for (var i = basisFunctions.length; i < values.length - 1; i++) {
        basisFunctions.add(BasisFunction(
            times[i - 1],
            times[i],
            times[i + 1],
            values[i]));
      }
      basisFunctions.add(BasisFunction(
          times[values.length - 2], times.last, null, values.last)); 
      // print(basisFunctions.length==values.length);
    }
  }

  Array interpolatePPG() {
    return interpolate(interpolation, millisInterp);
  }

  static Array interpolate(List<BasisFunction> basisFunctions, Array times) {
    var result = Array.empty();
    if (basisFunctions.isEmpty) return result;
    assert(times.first >= basisFunctions.first.midTime);
    assert(times.last <= basisFunctions.last.midTime);
    var indexInterp = 0;
    var indexTime = 0;
    while (indexTime < times.length) {
      if (times[indexTime] <= basisFunctions[indexInterp].finalTime) {
        result.add(basisFunctions[indexInterp].at(times[indexTime]) +
            basisFunctions[indexInterp + 1].at(times[indexTime]));
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

  static List findValleys(Array a, {double threshold}) {
    var N = a.length - 2;
    Array ix = Array.empty();
    Array ax = Array.empty();

    if (threshold != null) {
      for (int i = 1; i <= N; i++) {
        if (a[i - 1] >= a[i] && a[i] <= a[i + 1] && a[i] <= threshold) {
          ix.add(i.toDouble());
          ax.add(a[i]);
        }
      }
    } else {
      for (int i = 1; i <= N; i++) {
        if (a[i - 1] >= a[i] && a[i] <= a[i + 1]) {
          ix.add(i.toDouble());
          ax.add(a[i]);
        }
      }
    }
    return [ix, ax];
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

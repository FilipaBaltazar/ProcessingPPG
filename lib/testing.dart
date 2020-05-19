import 'package:smart_arrays_base/smart_arrays_base.dart';
import 'package:smart_arrays_numerics/smart_arrays_numerics.dart';
import 'package:smart_signal_processing/smart_signal_processing.dart';
import "dart:typed_data";
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';
import 'dart:math';
// import 'package:linalg/linalg.dart';

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
  static double eHbOxyRed = 298;

  /// Extinction coefficient of HbO at the green wavelength, in L / mmol / cm.
  static double eHbOxyGreen = 43016;

  /// Extinction coefficient of HbO at the blue wavelength, in L / mmol / cm.
  static double eHbOxyBlue = 34870.8;

  /// Extinction coefficient of Hb at the red wavelength, in L / mmol / cm.
  static double eHbRed = 2881.4;

  /// Extinction coefficient of Hb at the green wavelength, in L / mmol / cm.
  static double eHbGreen = 53412;

  /// Extinction coefficient of Hb at the blue wavelength, in L / mmol / cm.
  static double eHbBlue = 17025.6;

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
    var ratio = (eHbRed * sqrt(slopeBlue * peakBlue) -
            eHbBlue * sqrt(slopeRed * peakRed)) /
        (eHbOxyBlue * sqrt(slopeRed * peakRed) -
            eHbOxyRed * sqrt(slopeBlue * peakBlue));
    // var ratio = (eHbRed * (peakBlue) - eHbBlue * (peakRed)) /
    //     (eHbOxyBlue * (peakRed) - eHbOxyRed * (peakBlue));
    return 100 * (ratio / (ratio + 1));
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
        var peak2peak = valueHighPrev - valueLow;
        var deltaT = timeLow - timeHighPrev;
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

  /// [DateTime] values of each data point in the PPG signal.
  List<DateTime> dates = [];

  /// Timespan since [start] of each data point in the PPG signal, in milliseconds.
  Array millis = Array.empty();

  final List<BasisFunction> _interpolation = [];
  final Array _millisInterp = Array([0.0]);
  final Array _valuesInterp = Array.empty();
  final Array _valuesFiltered = Array.empty();
  final Array _valuesLog = Array.empty();
  final Map _peaks = {
    'indices': <int>[],
    'times': Array.empty(),
    'values': Array.empty(),
    'lastIndex': -1
  };
  final Map _valleys = {
    'indices': <int>[],
    'times': Array.empty(),
    'values': Array.empty(),
    'lastIndex': -1
  };
  double _pulseRate;
  final List<BasisFunction> _interpolationPeaks = [];
  final List<BasisFunction> _interpolationValleys = [];
  final Array _valuesMeanEnvelope = Array.empty();
  final Array _valuesProcessed = Array.empty();

  List<SensorValue> _sensorValuesRaw = [];
  List<SensorValue> _sensorValuesLog = [];
  List<SensorValue> _sensorValuesInterp = [];

  /// Sampling rate used to interpolate the PPG signal.
  int samplingRate;

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
    millis.add(time.difference(start).inMilliseconds.toDouble());
  }

  /// First value in [dates].
  DateTime get start {
    return dates.first;
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

  List<SensorValue> get sensorValuesRaw {
    var i = _sensorValuesRaw.length;
    while (i < valuesRaw.length) {
      _sensorValuesRaw.add(SensorValue(dates[i], valuesRaw[i]));
      i++;
    }
    return _sensorValuesRaw;
  }

  List<SensorValue> get sensorValuesLog {
    var i = _sensorValuesLog.length;
    while (i < valuesFiltered.length) {
      var date =
          start.add(Duration(microseconds: (millisInterp[i] * 1000).round()));
      _sensorValuesLog.add(SensorValue(date, valuesFiltered[i]));
      i++;
    }
    return _sensorValuesLog;
  }

  List<SensorValue> get sensorValuesInterp {
    var i = _sensorValuesInterp.length;
    while (i < valuesInterp.length) {
      var date =
          start.add(Duration(microseconds: (millisInterp[i] * 1000).round()));
      _sensorValuesInterp.add(SensorValue(date, valuesInterp[i]));
      i++;
    }
    return _sensorValuesInterp;
  }

  /// Returns the list of linear basis functions used to interpolate [valuesRaw].
  List<BasisFunction> get interpolation {
    if (_interpolation.length < valuesRaw.length) {
      _fillBasisFunctions(_interpolation, millis, valuesRaw);
    }
    return _interpolation;
  }

  /// Returns the length of [valuesInterp].
  int get lengthInterp {
    if (millis.isEmpty) {
      return 0;
    } else {
      return (millis.last * samplingRate / 1000).floor() + 1;
    }
  }

  /// Returns the timespan since [start] of each data point in [valuesInterp], in milliseconds.
  Array get millisInterp {
    while (_millisInterp.length < lengthInterp) {
      _millisInterp.add(_millisInterp.last + 1 / samplingRate * 1000);
    }
    return _millisInterp;
  }

  /// Returns the interpolated values obtained from [valuesRaw] at a sampling rate of [samplingRate].
  Array get valuesInterp {
    while (_valuesInterp.length < lengthInterp) {
      var times = Array(millisInterp.sublist(_valuesInterp.length));
      var valuesNew = interpolate(interpolation, times);
      for (var i = 0; i < valuesNew.length; i++) {
        _valuesInterp.add(valuesNew[i]);
      }
    }
    return _valuesInterp;
  }

  /// Returns a filtered version of [valuesInterp].
  Array get valuesFiltered {
    while (_valuesFiltered.length < valuesInterp.length) {
      // I should probably save these variables
      var frequencyNyquist = 0.5 * samplingRate;
      var bandPassWindow = Array([frequencyLow, frequencyHigh]);
      var normalBandPassWindow =
          bandPassWindow / Array([frequencyNyquist, frequencyNyquist]);
      var filterOrder = 30;
      var filterCoeffs =
          firwin(filterOrder, normalBandPassWindow, pass_zero: false);

      // get the end bit of valuesInterp, plus a tail with the size of the filter order
      var valuesNew = valuesInterp.getRangeArray(
          max(0, _valuesFiltered.length - filterOrder - 1),
          valuesInterp.length);
      // filter what we took from valuesInterp
      valuesNew = lfilter(filterCoeffs, Array([1.0]), valuesNew);
      // the number of new samples to add to the filtered signal
      var numNewSamples = valuesInterp.length - _valuesFiltered.length;
      for (var i = valuesNew.length - numNewSamples;
          i < valuesNew.length;
          i++) {
        _valuesFiltered.add(valuesNew[i]);
      }
    }
    return _valuesFiltered;
  }

  /// Returns the natural logarithm of [valuesFiltered].
  Array get valuesLog {
    while (_valuesLog.length < valuesFiltered.length) {
      // add the next term
      _valuesLog.add(log(valuesFiltered[_valuesLog.length]));
    }
    return _valuesLog;
  }

  /// Returns the indices, times, and values of the peaks of [valuesFiltered].
  ///
  /// A value is considered a peak if both
  /// the previous and next value are smaller or equal to it.
  Map get peaks {
    while (_peaks['lastIndex'] < valuesFiltered.length - 1) {
      // only check parts you haven't checked yet
      int firstIndex = max(0, _peaks['lastIndex'] - 1);
      var aux = findPeaks(
          valuesFiltered.getRangeArray(firstIndex, valuesFiltered.length));
      for (var i = 0; i < aux[0].length; i++) {
        var newIndex = aux[0][i].round() + firstIndex;
        _peaks['indices'].add(newIndex);
        _peaks['times'].add(millisInterp[newIndex]);
        _peaks['values'].add(aux[1][i]);
      }
      _peaks['lastIndex'] = valuesFiltered.length - 1;
    }
    return _peaks;
  }

  /// Returns the indices, times, and values of the valleys of [valuesFiltered].
  ///
  /// A value is considered a valley if both
  /// the previous and next value are greater or equal to it.
  Map get valleys {
    while (_valleys['lastIndex'] < valuesFiltered.length - 1) {
      // only check parts you haven't checked yet
      int firstIndex = max(0, _valleys['lastIndex'] - 1);
      var aux = findValleys(
          valuesFiltered.getRangeArray(firstIndex, valuesFiltered.length));
      for (var i = 0; i < aux[0].length; i++) {
        var newIndex = aux[0][i].round() + firstIndex;
        _valleys['indices'].add(newIndex);
        _valleys['times'].add(millisInterp[newIndex]);
        _valleys['values'].add(aux[1][i]);
      }
      _valleys['lastIndex'] = valuesFiltered.length - 1;
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
      var timeSpan = peaks['times'].last - peaks['times'].first;
      _pulseRate = ((peaks['times'].length - 1) / (timeSpan / 1000)) * 60;
      print('Peak count pulse rate');
      print(_pulseRate);

      var intervalsRR = arrayDiff(peaks['times']);
      var ratesRR = Array.empty();
      intervalsRR.forEach((i) => {ratesRR.add(1 / i * 1000 * 60)});
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

  List<BasisFunction> get interpolationPeaks {
    if (_interpolationPeaks.length < peaks['values'].length) {
      _fillBasisFunctions(_interpolationPeaks, peaks['times'], peaks['values']);
    }
    return _interpolationPeaks;
  }

  List<BasisFunction> get interpolationValleys {
    if (_interpolationValleys.length < valleys['values'].length) {
      _fillBasisFunctions(
          _interpolationValleys, valleys['times'], valleys['values']);
    }
    return _interpolationValleys;
  }

  /// Returns the first index of [valuesFiltered] that is covered by both envelopes.
  int get firstIndexEnvelopes {
    return max(peaks['indices'].first, valleys['indices'].first);
  }

  /// Returns the last index of [valuesFiltered] that is covered by both envelopes.
  int get lastIndexEnvelopes {
    return min(peaks['indices'].last, valleys['indices'].last);
  }

  /// Returns the length of [valuesMeanEnvelope], which is only defined
  /// in the region covered by both envelopes.
  int get lengthEnvelopes {
    return lastIndexEnvelopes - firstIndexEnvelopes + 1;
  }

  /// Returns the timepoints in [millisInterp] covered by both envelopes.
  Array get millisEnvelopes {
    return millisInterp.getRangeArray(
        firstIndexEnvelopes, lastIndexEnvelopes + 1);
  }

  /// Returns the mean envelope of [valuesFiltered].
  ///
  /// This array is only defined in the region covered by both the top and bottom envelopes,
  /// which is the region covered by both the peaks and valleys.
  Array get valuesMeanEnvelope {
    while (_valuesMeanEnvelope.length < lengthEnvelopes) {
      // lets only grab the times we don't yet have
      var times = Array(millisEnvelopes.sublist(_valuesMeanEnvelope.length));
      var valuesUpper = interpolate(interpolationPeaks, times);
      var valuesLower = interpolate(interpolationValleys, times);
      for (var i = 0; i < valuesUpper.length; i++) {
        _valuesMeanEnvelope.add(valuesUpper[i] - valuesLower[i]);
      }
    }
    return _valuesMeanEnvelope;
  }

  /// Returns the difference between [valuesFiltered] and [valuesMeanEnvelope].
  ///
  /// This array is only defined in the region covered by both the top and bottom envelopes,
  /// which is the region covered by both the peaks and valleys.
  Array get valuesProcessed {
    while (_valuesProcessed.length < lengthEnvelopes) {
      var newIndex = _valuesProcessed.length;
      // add the next entry
      _valuesProcessed.add(valuesFiltered[newIndex + firstIndexEnvelopes] -
          valuesMeanEnvelope[newIndex]);
    }
    return _valuesProcessed;
  }

  /// Fills `basisFunction` with [BasisFunction] objects
  /// until its length equals that of `values`.
  /// The resulting list serves as an interpolation function,
  /// where the images of the entries of `times` are the entries of `values`.
  static void _fillBasisFunctions(
      List<BasisFunction> basisFunctions, Array times, Array values) {
    if (values.isEmpty) {
      return;
    } else {
      // if basisFucntions is empty add the first (it's missing the startTime and the finalTime)
      if (basisFunctions.isEmpty) {
        basisFunctions.add(BasisFunction(null, times[0], null, values[0]));
      }

      // if basisFucntions is smaller than values, add until it's not
      if (basisFunctions.length < values.length) {
        // remove the last one because it's missing the finalTime
        basisFunctions.removeLast();

        // if length was 1, add back the first (do this seperately because startTime is null)
        if (basisFunctions.isEmpty) {
          basisFunctions
              .add(BasisFunction(null, times[0], times[1], values[0]));
        }

        // add functions until basisFunctions is one smaller than values
        while (basisFunctions.length < values.length - 1) {
          var i = basisFunctions.length;
          basisFunctions.add(
              BasisFunction(times[i - 1], times[i], times[i + 1], values[i]));
        }

        // add the last (it's missing the finalTime)
        var i = basisFunctions.length;
        basisFunctions
            .add(BasisFunction(times[i - 1], times[i], null, values[i]));
      }
    }
  }

  /// Returns the values given by the interpolation `basisFunctions` at the instants in `times`.
  static Array interpolate(List<BasisFunction> basisFunctions, Array times) {
    var result = Array.empty();

    if (basisFunctions.isEmpty) return result;

    assert(times.first >= basisFunctions.first.midTime,
        'The first time must be greater than the middle time of the first basis function.');
    assert(times.last <= basisFunctions.last.midTime,
        'The last time must be smaller than the middle time of the last basis function.');

    // find the index of the first basis function
    // whose middle time is smaller than the first entry of `times`
    var indexInterp = basisFunctions
        .indexWhere((function) => times.first >= function.midTime);
    // select the first time
    var indexTime = 0;

    while (indexTime < times.length) {
      // if the current time is smaller than
      // the final time of the current basis function,
      // then it is in that function's domain.
      // This is true because we have made sure that the first time is
      // greater than the middle time of the first basis function we selected,
      // and because the times and basisFunctions are sorted.
      // Thus, calculate the interpolated value and advance to the next time
      if (times[indexTime] <= basisFunctions[indexInterp].finalTime) {
        result.add(basisFunctions[indexInterp].at(times[indexTime]) +
            basisFunctions[indexInterp + 1].at(times[indexTime]));
        indexTime++;
      }
      // otherwise, the current time is not in the domain,
      // and we must increment the basis function
      else {
        indexInterp++;
      }
    }
    return result;
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

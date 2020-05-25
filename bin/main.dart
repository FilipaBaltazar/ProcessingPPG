import 'package:testing/testing.dart';
import 'package:scidart/numdart.dart';
import 'package:scidart/scidart.dart';
import 'dart:math';


class TestPLL {
static final Array _theta = Array.fixed(1000, initialValue: 0);
static final double phi = pi/12;
static final Array __p = Array([1, 2*cos(phi)-1]);

static Array get _p {
  while (__p.length < _theta.length) {
    var n = __p.length;
    __p.add(2*cos(_theta[n]+phi)*__p[n-1] - __p[n-2]);
  }
  return __p;
}

}


void main(List<String> arguments) {
  // var X = <double>[0,1,2,3,4,5];
  var times = [0.0, 0.079, 0.133, 0.187, 0.254, 0.319, 0.388, 0.441, 0.5, 0.564, 0.629, 1.012, 1.05, 1.105, 1.137, 1.183, 1.217, 1.248, 1.284, 1.317, 1.365, 1.39, 1.434, 1.466, 1.518, 1.574, 1.608, 1.638, 1.682, 1.712, 1.746, 1.786, 1.82, 1.878, 1.93, 1.968, 1.999, 2.033, 2.067, 2.105, 2.137, 2.18, 2.21, 2.247, 2.28, 2.316, 2.353, 2.385, 2.434, 2.464, 2.497, 2.53, 2.569, 2.614, 2.637, 2.68, 2.718, 2.746, 2.783, 2.813, 2.85, 2.885, 2.928, 2.96, 2.999, 3.033, 3.064, 3.103, 3.148, 3.178, 3.214, 3.246, 3.279, 3.315, 3.351, 3.384, 3.426, 3.46, 3.495, 3.528, 3.564, 3.601, 3.634, 3.675, 3.708, 3.744, 3.779, 3.815, 3.85, 3.891, 3.926, 3.959, 3.999, 4.029, 4.067, 4.1, 4.147, 4.172, 4.21, 4.242, 4.277, 4.312, 4.358, 4.392, 4.426, 4.458, 4.49, 4.528, 4.564, 4.599, 4.642, 4.674, 4.707, 4.742, 4.777, 4.815, 4.857, 4.891,
4.928, 4.959, 4.991, 5.028, 5.061, 5.106, 5.14, 5.173, 5.206, 5.242, 5.279, 5.314, 5.356, 5.39, 5.423, 5.456, 5.49, 5.527, 5.564, 5.598, 5.639, 5.672, 5.709, 5.743, 5.778, 5.814, 5.855, 5.889, 5.923, 5.957, 5.991, 6.027, 6.062, 6.096, 6.138, 6.173, 6.206, 6.242, 
6.276, 6.314, 6.355, 6.389, 6.423, 6.457, 6.489, 6.527, 6.563, 6.597, 6.637, 6.671, 6.706, 6.741, 6.777, 6.82, 6.857, 6.9, 6.923, 6.959, 6.998, 7.026, 7.062, 7.103, 7.881, 8.046, 8.12, 8.197, 
8.266, 8.359, 8.444, 8.513, 8.581, 8.656, 8.722, 8.785, 8.849, 8.907, 8.972, 9.04, 9.106, 9.179, 9.259, 9.28, 9.333, 9.395, 9.461, 9.527, 9.599, 9.673, 9.751, 9.814, 9.871, 9.925, 9.985, 10.037, 10.086, 10.152, 10.207, 10.268, 10.326, 10.388, 10.454, 10.519, 10.569, 10.622, 10.687, 10.705, 10.753, 10.808, 10.874, 10.931, 10.982, 11.038,
11.1, 11.155, 11.221, 11.241, 11.29, 11.354, 11.417, 11.477, 11.549, 11.611, 
11.671, 11.74, 11.813, 11.877, 11.945, 12.008, 12.074, 12.139, 12.215, 12.274, 12.336, 12.398, 12.454, 12.505, 12.569, 12.638, 12.7, 12.767, 12.834, 12.89, 12.949, 13.016, 13.067, 13.119, 13.181, 13.236, 13.298, 13.356, 13.418, 13.482, 13.533, 13.587, 13.65, 13.703, 13.766, 13.83, 13.886, 13.952, 14.016, 14.082, 14.136, 14.201, 14.278, 14.348, 14.41, 14.485, 14.571, 14.656, 14.718, 14.736, 14.794, 14.853, 14.916, 14.985, 15.053, 15.115, 15.182, 15.265, 15.343, 15.412, 15.479, 15.537, 15.601, 15.667, 15.737, 15.8, 15.865, 15.933, 15.95, 16.007, 16.067, 16.136, 16.199, 16.267, 16.333, 16.4, 16.472, 16.536, 16.601, 16.669, 16.734, 16.807, 16.867, 16.935, 17.007, 17.072, 17.135, 17.205, 17.268, 17.336, 17.4, 17.47, 17.536,
17.603, 17.671, 17.736, 17.802, 17.876, 17.938, 18.005, 18.073, 18.14, 
18.204, 18.268, 18.34, 18.403, 18.47, 18.54, 18.606, 18.673, 18.743, 18.81, 18.887, 18.957, 19.028, 19.1, 19.17, 19.237, 19.305, 19.372, 19.473, 19.539, 19.613, 19.69, 19.756, 19.824, 19.892, 
19.957, 20.025, 20.09, 20.157, 20.225, 20.298, 20.379, 20.457, 20.543, 20.621, 20.688, 20.759, 20.843, 20.925, 21.01, 21.088, 21.171, 21.254, 21.332, 21.43, 21.547, 21.627, 21.708, 21.803, 21.89, 
21.976, 22.083, 22.181, 22.284, 22.37, 22.457, 22.568, 22.657, 22.743, 22.827, 22.917, 22.994, 23.089, 23.167, 23.246, 23.33, 23.415, 23.493, 23.583, 23.669, 23.744, 23.829, 23.914, 23.995, 24.081, 24.166, 24.243, 24.329, 24.415, 24.497, 24.579, 24.666, 24.75, 24.829, 24.916, 25.001, 25.101, 25.188, 25.28, 25.364, 25.457, 25.551, 25.652, 25.759, 25.861, 25.961, 26.055, 26.154, 26.253, 26.358, 26.456, 26.555,
26.657, 26.752, 26.854, 26.957, 27.055, 27.155, 27.261, 27.357, 27.46, 27.573, 27.683, 27.789, 27.898, 27.994, 28.092, 28.196, 28.292, 28.394, 28.493, 28.591, 28.699, 28.812, 28.913, 29.024, 29.135, 29.247, 29.265, 29.363, 29.468, 29.575, 29.682, 29.801, 29.913, 30.018, 30.129, 30.245, 30.371, 30.499, 30.632, 30.776, 30.903, 31.043, 31.189, 31.329, 31.468, 31.61, 31.748, 31.766, 31.896, 32.04, 32.164, 32.293, 32.416, 32.547, 32.676, 32.81, 32.873, 33.015];

  var pipa_heart = [251.67619791666667, 253.5387109375, 253.50421875, 237.00524739583332, 235.77805989583334, 96.7458984375, 85.16903645833334, 90.21041666666666, 90.5400390625, 114.19184895833334, 134.82263020833332, 141.74333333333334, 148.23321614583332, 144.63485677083332, 143.74013020833334, 142.7415625, 144.17709635416668, 145.9291015625, 148.3608984375, 148.18438802083332, 150.27565104166666, 150.36287760416667, 152.01627604166666, 153.23596354166668, 156.66645833333334, 159.666171875, 160.3045703125, 163.06079427083333, 163.96856770833332, 165.10572916666666, 165.9142578125, 162.51143229166667, 163.11623697916667, 162.2674609375, 161.92283854166666, 162.35891927083333, 162.13311197916667, 161.49108072916667, 163.55467447916666, 158.707578125, 157.57364583333333, 156.93322916666668, 160.02282552083332, 160.41546875, 161.9524609375, 164.18666666666667, 164.91729166666667, 164.94888020833332, 165.70244791666667, 165.78006510416665, 165.54794270833332, 165.60337239583333, 165.83884114583333, 165.917421875, 167.02677083333333, 167.45489583333332, 167.60692708333335, 167.65998697916666, 167.84029947916667, 168.03444010416666, 168.07463541666667, 167.63451822916667, 165.67666666666668, 163.2169921875, 162.00360677083333, 162.00541666666666, 162.60825520833333, 163.15979166666668, 163.59557291666667, 163.82471354166665, 163.7102734375, 163.4027734375, 163.18557291666667, 163.35467447916668, 163.65850260416667, 163.8984765625, 163.69638020833332, 163.27864583333334, 162.85092447916668, 162.57061197916667, 162.07446614583333, 161.88666666666666, 161.5881640625, 161.46404947916668, 161.08703125,161.08122395833334, 160.83782552083332, 160.91953125, 160.57052083333335, 159.09454427083332, 156.45036458333334, 155.10165364583332, 154.47583333333333, 154.467421875, 154.45130208333333, 154.56653645833333, 154.80611979166667, 154.89553385416667, 154.88403645833333, 154.8966015625, 155.07173177083334, 155.4070703125, 155.85295572916667, 156.3980078125, 156.82802083333334, 157.26002604166666, 157.6766796875, 158.096875, 158.5266796875, 158.79342447916667, 159.2661328125, 159.582421875, 160.14973958333334, 160.45770833333333, 160.927890625, 160.68981770833332, 159.186484375, 157.24833333333333, 159.167890625, 159.81122395833333, 160.07635416666668, 163.1566796875, 164.285625, 164.4671484375, 165.22932291666666, 165.42697916666665, 165.470078125, 165.93118489583333, 166.3442578125, 166.70591145833333, 167.10069010416666, 167.5171484375, 167.93783854166668, 168.2347265625, 168.38315104166668, 168.64209635416665, 168.79783854166666, 169.01196614583333, 169.18390625, 169.40381510416665, 169.05670572916668, 167.3158984375, 165.70725260416665, 165.0251171875, 165.10106770833335, 165.25946614583333, 165.55080729166667, 165.72709635416666, 166.04403645833332, 166.07522135416667, 166.0216015625, 165.95786458333333, 166.1675390625, 166.4678515625, 166.89248697916668, 167.290078125, 167.708046875, 167.9877734375, 168.19713541666667, 168.40115885416665, 168.4759765625, 168.54557291666666, 168.64376302083335, 168.94231770833332, 169.08458333333334, 169.3946484375, 169.65377604166667, 169.73609375, 168.40739583333334, 166.3276953125, 164.98311197916667, 164.49389322916667, 164.60837239583333, 164.62067708333333, 164.75052083333333, 164.93256510416666, 165.22401041666666, 165.19709635416666, 165.14877604166668, 165.15095052083333, 165.44815104166668, 167.88997395833334, 161.32588541666667, 163.60930989583332, 162.11694010416667, 159.90494791666666, 158.27651041666667, 163.17986979166668, 163.51484375, 163.69786458333334, 164.09192708333333, 164.28135416666666, 164.73170572916666, 165.14880208333332, 165.60973958333332, 163.63513020833332, 161.44703125, 160.18826822916665, 160.73928385416667, 161.43151041666667, 161.52509114583333, 161.34692708333333, 161.1325, 161.86397135416667, 162.67770833333333, 163.2482421875, 163.63372395833332, 163.9675390625, 164.41274739583332, 164.94790364583332, 165.19765625, 165.73776041666667, 166.0019140625, 165.16165364583333, 160.72477864583334, 160.00115885416668, 160.2980078125, 160.52361979166668, 161.1758203125, 160.8650390625, 160.90143229166668, 161.24653645833334, 162.133046875, 162.4448046875, 162.69959635416666, 162.935078125, 163.57411458333334, 164.07506510416667, 164.42571614583332, 164.94403645833333, 165.13532552083333, 165.15671875, 165.44890625, 165.73713541666666, 166.04272135416667, 166.183359375, 162.49783854166665, 160.48493489583333, 161.17536458333333, 162.08631510416666, 162.3694140625, 161.88265625, 162.08173177083333, 163.03298177083335, 163.70768229166666, 164.00786458333334, 164.240234375, 164.6171484375, 164.9882421875, 165.48811197916666, 166.02897135416666, 163.7965625, 160.60635416666668, 160.67930989583334, 161.08401041666667, 161.98911458333333, 162.21548177083332, 161.751015625, 162.19311197916667,
 162.68471354166667, 163.43092447916666, 163.8266015625, 164.12364583333334, 164.28579427083332, 164.60576822916667, 164.80803385416667, 165.31407552083334, 165.2976953125, 163.38182291666666, 159.90236979166667, 159.9596484375, 160.22751302083333, 160.89408854166666, 161.10234375, 160.89684895833332, 161.09063802083332, 162.039765625, 162.48513020833335, 163.12326822916665, 163.564296875, 163.8192578125, 164.09735677083333, 164.5480859375, 165.018125, 165.59516927083334, 164.604296875, 160.178984375, 159.63505208333333, 160.45109375, 160.82108072916665, 161.0598046875, 160.92770833333333, 161.09796875, 161.92416666666668, 162.84588541666668, 163.33166666666668, 163.84700520833334, 164.0155078125, 164.66510416666668, 165.03873697916666, 165.48216145833334, 165.96069010416667, 166.15794270833334,
 165.78138020833333, 161.1200390625, 160.30549479166666, 161.3301953125, 162.1960546875, 162.21828125, 161.92973958333334, 161.66509114583334, 162.07252604166666, 162.9844140625, 163.53438802083335, 163.9275, 164.18184895833335, 164.60768229166666, 165.08529947916668, 165.31916666666666, 165.03682291666667, 160.8759375, 160.42463541666666, 161.71981770833332, 162.60459635416666, 162.11809895833332, 162.323046875, 162.80325520833333, 163.61397135416667, 164.14567708333334, 164.33283854166666, 164.6050390625, 165.10140625, 165.34600260416667, 163.885703125, 159.90462239583334, 160.02259114583333, 160.9402734375, 161.65729166666668, 161.26838541666666, 161.67489583333332, 162.16149739583332, 163.01307291666666, 163.6671875, 164.08322916666665, 164.37498697916666, 164.8630078125,
165.35700520833333, 164.44149739583332, 159.99447916666668, 159.47908854166667, 160.17063802083334, 160.965078125, 160.81997395833332, 161.05319010416667, 162.02665364583333, 162.82954427083334, 163.35799479166667, 163.71673177083332, 164.0474609375, 164.45157552083333, 164.9780078125, 163.16877604166666, 160.59915364583333, 159.52846354166667, 160.33010416666667, 160.98548177083333, 160.88634114583334, 160.77959635416667, 161.64817708333334, 162.48662760416667, 163.13033854166667, 163.64779947916668, 163.9045703125, 164.28772135416668, 164.75309895833334, 165.19384114583335, 161.0679296875, 159.48514322916665, 160.52333333333334, 161.63713541666667, 161.4034375, 161.31411458333332, 162.71490885416668, 163.36231770833334, 163.86240885416666, 164.35619791666667, 164.86751302083334, 165.09904947916667,
159.40393229166668, 160.15920572916667, 161.13873697916668, 160.810078125, 161.30368489583333, 162.711171875, 163.52143229166666, 163.8848046875, 164.42006510416667, 164.48674479166667, 160.04274739583335, 158.89608072916667, 159.7454296875, 159.83307291666668, 160.10450520833334, 161.01709635416665, 161.7547265625, 162.596328125, 163.01580729166668, 163.54669270833332, 164.3830078125, 162.11819010416667, 158.47592447916668, 158.66572916666667, 159.36080729166667, 159.45815104166667, 159.63244791666668, 160.5915234375, 161.3709375, 162.189296875, 162.6253125, 163.07846354166668, 163.84958333333333, 163.77071614583335, 159.6079296875, 159.02510416666667, 159.92576822916666, 160.21373697916667, 159.92544270833332, 160.76958333333334, 161.8551171875, 162.27720052083333, 162.8183203125, 
163.33984375, 161.87809895833334, 158.6859765625, 158.79315104166668, 
159.5956640625, 159.29157552083333, 160.33537760416667, 160.96911458333332, 161.69438802083334, 
162.60549479166667, 162.24514322916667, 158.51430989583332, 157.61875, 158.18869791666665, 157.98401041666668, 158.90018229166665, 160.12709635416667, 160.68817708333333, 161.5430859375, 162.43842447916666, 163.38890625, 158.61834635416668, 157.4928125, 157.73390625, 157.89319010416668, 
158.789453125, 159.59013020833333, 160.67411458333333, 161.67303385416668, 162.5527734375, 163.30348958333335, 162.55888020833333, 158.82330729166668, 158.54434895833333, 159.15680989583333, 158.746640625, 159.71583333333334, 160.69868489583334, 161.39563802083333, 162.2615625, 162.629453125, 161.9409765625, 157.79752604166666, 158.07891927083332,
158.32372395833335, 158.39514322916668, 159.50766927083333, 160.59111979166667, 161.50240885416667, 162.415703125, 162.72825520833334, 157.5253125, 157.81717447916665, 157.8810546875, 159.1750390625, 160.56055989583334, 161.79977864583333, 163.05657552083332, 161.21270833333332, 158.18828125, 158.99346354166667, 159.1014453125, 159.109765625, 160.65110677083334, 161.35790364583335, 162.4768359375, 163.685390625, 159.63348958333333, 160.631171875, 163.408671875, 162.3865234375];

  var last = DateTime.now();
  var start= last;

  PPG ppg = PPG(25);


  for (var i = 0; i < pipa_heart.length; i++) {
    last = start.add(Duration(milliseconds:(times.elementAt(i)*1000).round()));
  
    ppg.add(last, pipa_heart.elementAt(i));  

    print(ppg.valuesProcessed.length);

    // print(last.difference(start).inMilliseconds);  
  }

  ppg.pulseRate;

  //print(findPeaks(Array([1,2,3]))[0]);

  //PPG test = PPG(25);

  // var timesNew = [0, 1, 2, 3, 4, 5, 6, 7];
  // var valuesNew = <double>[3.5, 3, 2.5, 2, 1.5, 1, .5];
  // last = DateTime.now();
  // for (var i = 0; i <  valuesNew.length; i++) {
  //   last = start.add(Duration(milliseconds:(timesNew.elementAt(i)*1000).round()));
  //   test.add(last, valuesNew.elementAt(i)); 
  // }


  // print(Array(List.generate(ppg.millisInterp.length, (index) => ppg.window(ppg.millisInterp[index]/1000))));
  // print(ppg.reflectionPoint);
  // print(ppg.conservedIndex(ppg.millisInterp.length-1)/ppg.samplingRate);

  //print(ppg.pulseRate);
  //print(ppg.projections);
  // print(ppg.valuesProcessed);

  //print(ppg.projection(60-40));

  //print(test.pulseRate.toString());
  //print(ppg.projections);
  // print(test.valuesProcessed);

  // print(ppg.valuesInterp.length);

  // print(ppg.valuesInterp.length);

  //print(TestPLL._p);

}

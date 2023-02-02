#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>

#include <Garfield/AvalancheMC.hh>
#include <Garfield/ComponentAnalyticField.hh>
#include <Garfield/MediumDiamond.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/Shaper.hh>
#include <Garfield/TrackHeed.hh>
#include <Garfield/TrackTrim.hh>
#include <Garfield/ViewCell.hh>
#include <Garfield/ViewField.hh>
#include <Garfield/ViewSignal.hh>

int main(int argc, char *argv[]) {
  /**
   * @brief This routine simulates the sCVD signal
   *
   */
  using namespace Garfield;

  TApplication app("app", &argc, argv);

  MediumDiamond medium;
  medium.SetTemperature(300.0);

  // Define the cell layout.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&medium);
  cmp.AddPlaneY(0.0, 0., "p1");
  cmp.AddPlaneY(0.05, 400., "p2");

  // Add a readout strip along z.
  cmp.AddPixelOnPlaneY(0, -0.2, 0.2, -0.2, 0.2, "pixel");
  cmp.AddReadout("pixel");
  cmp.PrintCell();

  Shaper shaper(1, 80., 1, "unipolar");

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "pixel");
  sensor.SetTransferFunction(shaper);

  // Plot the weighting potential.
  const double xmin = -0.01;
  const double xmax = 0.01;
  const double ymin = 0.0;
  const double ymax = 0.05;
  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetSensor(&sensor);
  fieldView.SetArea(xmin, ymin, xmax, ymax);
  // fieldView.PlotContourWeightingField("pixel", "v");
  fieldView.Plot("v", "COLZ0");

  // Timming information
  const unsigned int nTimeBins = 5000;
  const double tmin = -100.;
  const double tmax = 400.;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Track TRIM
  TrackTrim track;
  track.SetSensor(&sensor);
  const std::string filename = "data/EXYZ.txt";
  const unsigned int nIons = 1;
  const unsigned int nSkip = 0;
  if (!track.ReadFile(filename, nIons, nSkip)) {
    std::cerr << "Reading TRIM EXYZ file failed.\n";
    return 1;
  }
  track.Print();

  // Simulate electron/hole drift lines using MC integration.
  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.EnableSignalCalculation();
  drift.EnableAttachmentMap();

  ViewDrift vDrift;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    vDrift.SetArea(xmin, ymin, xmax, ymax);
    track.EnablePlotting(&vDrift);
    drift.EnablePlotting(&vDrift);
  }

  double x0 = 0 * 1.5, y0 = 0.05, z0 = 0., t0 = 0.;
  double dx = 0., dy = -1., dz = 0.;
  double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
  int ne = 0;
  // Retrieve the clusters along the track.
  for (unsigned int i = 0; i < nIons; ++i) {
    track.NewTrack(x0, y0, z0, 0, dx, dy, dz);
    while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
      drift.SetElectronSignalScalingFactor(ne);
      drift.DriftElectron(xc, yc, zc, tc);
      drift.SetHoleSignalScalingFactor(ne);
      drift.DriftHole(xc, yc, zc, tc);
    }
  }

  // sensor.ConvoluteSignals();

  if (plotDrift) {
    TCanvas *cd = new TCanvas("cd", "", 600, 600);
    vDrift.SetPlane(0, 0, -1, 0, 0, 0);
    vDrift.SetArea(xmin, ymin, xmax, ymax);
    vDrift.SetCanvas(cd);
    constexpr bool twod = true;
    vDrift.Plot(twod);
  }

  constexpr bool plotSignal = true;
  ViewSignal signalView;
  signalView.SetSensor(&sensor);
  constexpr bool plotTotalSignal = true;
  constexpr bool plotElectronSignal = true;
  constexpr bool plotHoleSignal = true;
  signalView.PlotSignal("pixel", plotTotalSignal, plotElectronSignal,
                        plotHoleSignal, true);

  gSystem->ProcessEvents();
  app.Run(true);

  return 0;
}

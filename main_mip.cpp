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
  /* code */
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

  Shaper shaper(1, 80., 12.5, "unipolar");

  // Declare sensor
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

  // Set up Heed.
  TrackHeed track;
  track.SetSensor(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("muon");
  track.SetMomentum(100.e9);

  const unsigned int nMIP = 1;

  // Simulate electron/hole drift lines using MC integration.
  AvalancheMC drift;
  drift.SetSensor(&sensor);
  // Use steps of 1 micron.
  // drift.SetDistanceSteps(1.e-4);
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
  for (unsigned int i = 0; i < nMIP; ++i) {
    track.NewTrack(x0, y0, z0, 0, dx, dy, dz);

    while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
      //   // Loop over the electrons in the cluster.
      for (int j = 0; j < ne; ++j) {
        double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
        double dxe = 0., dye = 0., dze = 0.;
        track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
        // Simulate the electron and hole drift lines.
        drift.DriftElectron(xe, ye, ze, te);
        drift.DriftHole(xe, ye, ze, te);
      }
    }
  }

  //sensor.ConvoluteSignals();
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

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

void convert_to_root(const char *infile = "tracking.txt",
                     const char *outfile = "tracking_1GeV_Xe.root") {

  std::ifstream fin(infile);
  if (!fin.is_open()) {
    std::cerr << "Cannot open " << infile << std::endl;
    return;
  }

  TFile *fout = TFile::Open(outfile, "RECREATE");

  double x, y, z, Px, Py, Pz, t, Weight;
  int PDGid, EventID, TrackID, ParentID;

  std::map<std::string, TTree *> trees;

  const std::string detNames[] = {"DetIn1",  "DetIn2",  "DetIn3",  "DetIn4",
                                  "DetIn5",  "DetOut1", "DetOut2", "DetOut3",
                                  "DetOut4", "DetOut5"};

  // =========================
  // Detector centers (mm)
  // =========================
  struct DetPos {
    std::string name;
    double x;
    double z;
  };

  std::vector<DetPos> detPos = {
      {"DetIn1", 0.0, 200.0},      {"DetIn2", 0.0, 500.0},
      {"DetIn3", 0.0, 800.0},      {"DetIn4", 0.0, 1100.0},
      {"DetIn5", 0.0, 1400.0},

      {"DetOut1", 483.0, 2283.0},  {"DetOut2", 786.0, 2566.0},
      {"DetOut3", 1089.0, 2849.0}, {"DetOut4", 1392.0, 3132.0},
      {"DetOut5", 1695.0, 3415.0}};

  // =========================
  // Create trees
  // =========================
  for (auto &name : detNames) {
    TTree *tr = new TTree(name.c_str(), name.c_str());

    tr->Branch("x", &x);
    tr->Branch("y", &y);
    tr->Branch("z", &z);
    tr->Branch("Px", &Px);
    tr->Branch("Py", &Py);
    tr->Branch("Pz", &Pz);
    tr->Branch("t", &t);
    tr->Branch("PDGid", &PDGid);
    tr->Branch("EventID", &EventID);
    tr->Branch("TrackID", &TrackID);
    tr->Branch("ParentID", &ParentID);
    tr->Branch("Weight", &Weight);

    trees[name] = tr;
  }

  // =========================
  // Read ASCII
  // =========================
  std::string line;

  while (getline(fin, line)) {

    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream ss(line);

    if (!(ss >> x >> y >> z >> Px >> Py >> Pz >> t >> PDGid >> EventID >>
          TrackID >> ParentID >> Weight))
      continue;

    // =========================
    // Find nearest detector
    // =========================
    double dmin = std::numeric_limits<double>::max();
    std::string det;

    for (const auto &d : detPos) {

      double dx = x - d.x;
      double dz = z - d.z;

      double dist2 = dx * dx + dz * dz;

      if (dist2 < dmin) {
        dmin = dist2;
        det = d.name;
      }
    }

    trees[det]->Fill();
  }

  fout->Write();
  fout->Close();

  std::cout << "Convert: Done" << std::endl;
}

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TVectorD.h>

using ll = long long;

// ============================================================
// STRUCTURE MEASUREMENT FROM ROOT
// ============================================================
const int nDetIn = 5;
const int nDetOut = 5;
const double h_drift_in = 300.0; // mm

// const double h_B = 720.0;         // mm, arc length DetIn5→DetOut1
const double h_B = 1000.6;

// const double h_drift_out = 400.0; // mm
const double h_drift_out = 415;

// const double p = 230.0;           // MeV
const double p = 1000.6; // MeV

// const double B = 1.5;             // T
const double B = 6.54; // T

const double rho = (p / 1000.0) / (0.3 * B) * 1000.0;

const std::string detInName[nDetIn] = {"DetIn1", "DetIn2", "DetIn3", "DetIn4",
                                       "DetIn5"};
const std::string detOutName[nDetOut] = {"DetOut1", "DetOut2", "DetOut3",
                                         "DetOut4", "DetOut5"};

// ============================================================
// RK4 PREDICT
// State vector 4D: (y, ty, x, tx)
// state: (y, ty), dy/ds = ty, dty/ds = f2
// rho=0 => drift (f2=0), rho>0 => B field (f2=-(-1/rho))
// ============================================================
TVectorD RK4_predict(const TVectorD &x, double h, double rho_param = 0.0) {
  double f_ty = 0.0;                                       // ty linear tracking
  double f_tx = (rho_param > 0.0) ? 1.0 / rho_param : 0.0; // check -(-...)

  // k1
  double k1_y = x[1];
  double k1_ty = f_ty;
  double k1_x = x[3];
  double k1_tx = f_tx;

  // k2
  double k2_y = x[1] + h / 2.0 * k1_ty;
  double k2_ty = f_ty;
  double k2_x = x[3] + h / 2.0 * k1_tx;
  double k2_tx = f_tx;

  // k3
  double k3_y = x[1] + h / 2.0 * k2_ty;
  double k3_ty = f_ty;
  double k3_x = x[3] + h / 2.0 * k2_tx;
  double k3_tx = f_tx;

  // k4
  double k4_y = x[1] + h * k3_ty;
  double k4_ty = f_ty;
  double k4_x = x[3] + h * k3_tx;
  double k4_tx = f_tx;

  TVectorD xn(4);                                                   // 4D
  xn[0] = x[0] + h / 6.0 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);     // y
  xn[1] = x[1] + h / 6.0 * (k1_ty + 2 * k2_ty + 2 * k3_ty + k4_ty); // ty
  xn[2] = x[2] + h / 6.0 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);     // x
  xn[3] = x[3] + h / 6.0 * (k1_tx + 2 * k2_tx + 2 * k3_tx + k4_tx); // tx

  return xn;
}

// ============================================================
// Q MATRIX (4 X 4)
// ============================================================
TMatrixD make_Q(double h, double sigma_ms) {
  double s2 = sigma_ms * sigma_ms;
  TMatrixD Q(4, 4);
  Q.Zero();
  // y block
  Q(0, 0) = h * h * s2;
  Q(0, 1) = h * s2;
  Q(1, 0) = h * s2;
  Q(1, 1) = s2;
  // x block
  Q(2, 2) = h * h * s2;
  Q(2, 3) = h * s2;
  Q(3, 2) = h * s2;
  Q(3, 3) = s2;

  return Q;
}

// ============================================================
// STRUCTURE KF
// ============================================================
struct StepResult {
  double res_y;
  double res_x;
  double S_yy;
  double S_xx;
};

StepResult KF_step(TVectorD &x, TMatrixD &Cov, double h, double rho_param,
                   double y_meas_val, double x_meas_val, double sigma_det,
                   const TMatrixD &Q) {
  // PREDICT
  TVectorD x_pred = RK4_predict(x, h, rho_param);

  TMatrixD F(4, 4);
  F.Zero();
  F(0, 0) = 1;
  F(0, 1) = h;
  F(1, 1) = 1;
  F(2, 2) = 1;
  F(2, 3) = h;
  F(3, 3) = 1;

  TMatrixD Cov_pred = F * Cov * TMatrixD(TMatrixD::kTransposed, F) + Q;

  // UPDATE
  TMatrixD H(2, 4);
  H.Zero();
  H(0, 0) = 1; // y
  H(1, 2) = 1; // x

  TMatrixD R(2, 2);
  R.Zero();
  R(0, 0) = sigma_det * sigma_det;
  R(1, 1) = sigma_det * sigma_det;

  TVectorD z(2);
  z[0] = y_meas_val;
  z[1] = x_meas_val;

  TVectorD res_vec = z - H * x_pred;

  TMatrixD S = H * Cov_pred * TMatrixD(TMatrixD::kTransposed, H) + R;
  TMatrixD S_inv = S;
  S_inv.Invert();

  TMatrixD K = Cov_pred * TMatrixD(TMatrixD::kTransposed, H) * S_inv;
  x = x_pred + K * res_vec;

  TMatrixD I(4, 4);
  I.UnitMatrix();
  TMatrixD ImKH = I - K * H;
  Cov = ImKH * Cov_pred * TMatrixD(TMatrixD::kTransposed, ImKH) +
        K * R * TMatrixD(TMatrixD::kTransposed, K);

  return {res_vec[0], res_vec[1], S(0, 0), S(1, 1)};
}

// ============================================================
// MAIN
// ============================================================
void kf_rk4() {

  TFile *file = new TFile("tracking_1GeV_Xe.root");
  TRandom3 rng(42);
  double sigma_det = 0.01;
  double sigma_ms = 0.016; // = sigma_ms_Xe_1Gauss

  // Finding data base
  std::vector<std::vector<double>> y_meas_in(nDetIn), y_true_in(nDetIn);
  std::vector<std::vector<double>> x_meas_in(nDetIn), x_true_in(nDetIn);
  std::vector<std::vector<double>> y_meas_out(nDetOut), y_true_out(nDetOut);
  std::vector<std::vector<double>> x_meas_out(nDetOut), x_true_out(nDetOut);

  for (int i = 0; i < nDetIn; i++) {
    TTree *t = (TTree *)file->Get(detInName[i].c_str());
    double y, x_pos;
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("x", &x_pos);
    for (ll e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      y_true_in[i].push_back(y);
      x_true_in[i].push_back(x_pos);
      y_meas_in[i].push_back(y + rng.Gaus(0, sigma_det));
      x_meas_in[i].push_back(x_pos + rng.Gaus(0, sigma_det));
    }
  }
  for (int i = 0; i < nDetOut; i++) {
    TTree *t = (TTree *)file->Get(detOutName[i].c_str());
    double y, x_pos;
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("x", &x_pos);
    for (ll e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      y_true_out[i].push_back(y);
      x_true_out[i].push_back(x_pos);
      y_meas_out[i].push_back(y + rng.Gaus(0, sigma_det));     // + ms
      x_meas_out[i].push_back(x_pos + rng.Gaus(0, sigma_det)); // + ms
    }
  }

  ll nEvents = y_meas_in[0].size();

  const int ndof_y = nDetIn + nDetOut;
  const int ndof_x = nDetIn + nDetOut;

  // ============================================================
  // HISTOGRAMS
  // ============================================================
  TH1D *h_pull_y = new TH1D("h_pull_y", "Pull Y;pull_y;counts", 100, -5, 5);
  TH1D *h_pull_x = new TH1D("h_pull_x", "Pull X;pull_x;counts", 100, -5, 5);
  TH1D *h_chi2_y = new TH1D("h_chi2_y", "Chi2 Y;#chi^{2};counts", 300, 0, 50);
  TH1D *h_chi2_x = new TH1D("h_chi2_x", "Chi2 X;#chi^{2};counts", 300, 0, 50);
  TH1D *h_pval =
      new TH1D("h_pval", "P Value Total;p_val;counts", 300, -0.05, 1.05);

  // HISTOGRAMS - RESULT
  TH1D *h_momentum_results = new TH1D(
      "h_momentum_results", "Predict Momentum;p [GeV/c];counts", 300, 0.5, 1.5);
  TH1D *h_dp =
      new TH1D("h_dp", "Momentum Resolution;#Deltap/p;counts", 300, -0.5, 0.5);

  double h_x_B = x_true_out[0][0] - x_true_in[4][0];   // DetIn5 -> DetOut1
  double dx_out = x_true_out[1][0] - x_true_out[0][0]; // DetOut1 -> DetOut2

  // ============================================================
  // KF EVENT LOOP
  // ============================================================
  TVectorD x(4);
  TMatrixD Cov(4, 4);

  // Q Matrix
  TMatrixD Q_in = make_Q(h_drift_in, sigma_ms);
  TMatrixD Q_B = make_Q(h_B, sigma_ms);
  TMatrixD Q_out = make_Q(h_drift_out, sigma_ms);

  for (ll evt = 0; evt < nEvents; evt++) {
    // Init for KF
    x[0] = y_meas_in[0][evt];
    x[1] = (y_meas_in[1][evt] - y_meas_in[0][evt]) / h_drift_in;
    x[2] = x_meas_in[0][evt];
    x[3] = (x_meas_in[1][evt] - x_meas_in[0][evt]) / h_drift_in;

    Cov.Zero();
    double s2 = sigma_det * sigma_det;
    double h2 = h_drift_in * h_drift_in;
    // Covariance Y
    Cov(0, 0) = s2;
    Cov(0, 1) = -s2 / h_drift_in;
    Cov(1, 0) = -s2 / h_drift_in;
    Cov(1, 1) = 2 * s2 / h2;
    // Covariance X
    Cov(2, 2) = s2;
    Cov(2, 3) = -s2 / h_drift_in;
    Cov(3, 2) = -s2 / h_drift_in;
    Cov(3, 3) = 2 * s2 / h2;
    double chi2_x = 0.0;
    double chi2_y = 0.0;
    {
      // Pull det1
      double pull_y = (x[0] - y_true_in[0][evt]) / sqrt(Cov(0, 0));
      double pull_x = (x[2] - x_true_in[0][evt]) / std::sqrt(Cov(2, 2));
      h_pull_y->Fill(pull_y);
      h_pull_x->Fill(pull_x);

      // chi2
      chi2_y += std::pow((x[0] - y_true_in[0][evt]) / sqrt(Cov(0, 0)), 2);
      chi2_x += std::pow((x[2] - x_true_in[0][evt]) / std::sqrt(Cov(2, 2)), 2);
    }

    // ==============================
    // DetIn2 => DetIn5
    // ==============================
    for (int i = 1; i < nDetIn; i++) {
      auto [ry, rx, Syy, Sxx] =
          KF_step(x, Cov, h_drift_in, 0.0, y_meas_in[i][evt], x_meas_in[i][evt],
                  sigma_det, Q_in);

      // Tracking Y
      double pull_y = (x[0] - y_true_in[i][evt]) / std::sqrt(Cov(0, 0));
      double pull_x = (x[2] - x_true_in[i][evt]) / std::sqrt(Cov(2, 2));

      h_pull_y->Fill(pull_y);
      h_pull_x->Fill(pull_x);
      chi2_y += std::pow((x[0] - y_true_in[i][evt]), 2) / Cov(0, 0);
      chi2_x += std::pow((x[2] - x_true_in[i][evt]), 2) / Cov(2, 2);
    }

    // CHECK ERROR

    double tx_in = x[3];
    // ==============================
    // DetIn5 → DetOut1
    // ==============================
    {
      auto [ry, rx, Syy, Sxx] = KF_step(x, Cov, h_B, rho, y_meas_out[0][evt],
                                        x_meas_out[0][evt], sigma_det, Q_B);

      // pull
      double pull_y = (x[0] - y_true_out[0][evt]) / std::sqrt(Cov(0, 0));
      double pull_x = (x[2] - x_true_out[0][evt]) / std::sqrt(Cov(2, 2));
      h_pull_y->Fill(pull_y);
      h_pull_x->Fill(pull_x);

      // chi2
      chi2_y += std::pow((x[0] - y_true_out[0][evt]), 2) / Cov(0, 0);
      chi2_x += std::pow((x[2] - x_true_out[0][evt]), 2) / Cov(2, 2);
    }
    double tx_out = x[3];
    double theta_in = atan(tx_in);
    double theta_out = atan(tx_out);
    double delta_theta = fabs(theta_out - theta_in); // saggita
    double L_mag =
        (400 / std::cos(30 * TMath::Pi() / 180)) / 1000.0; // Length of magnetic
    double momen_predict = 0.3 * B * L_mag / delta_theta;  // P (GeV)
    h_momentum_results->Fill(momen_predict);               // p_true = 1 GeV/c
    h_dp->Fill((momen_predict - 1.0) / 1.0);

    // CHECK ERROR

    // ==============================
    // DetOut2 → DetOut5
    // ==============================
    for (int i = 1; i < nDetOut; i++) {
      auto [ry, rx, Syy, Sxx] =
          KF_step(x, Cov, h_drift_out, 0.0, y_meas_out[i][evt],
                  x_meas_out[i][evt], sigma_det, Q_out);

      // pull
      double pull_y = (x[0] - y_true_out[i][evt]) / std::sqrt(Cov(0, 0));
      double pull_x = (x[2] - x_true_out[i][evt]) / std::sqrt(Cov(2, 2));
      h_pull_y->Fill(pull_y);
      h_pull_x->Fill(pull_x);

      // chi2
      chi2_y += std::pow((x[0] - y_true_out[i][evt]), 2) / Cov(0, 0);
      chi2_x += std::pow((x[2] - x_true_out[i][evt]), 2) / Cov(2, 2);
    }

    h_chi2_y->Fill(chi2_y);
    h_chi2_x->Fill(chi2_x);
    h_pval->Fill(TMath::Prob(chi2_y + chi2_x,
                             ndof_x + ndof_y)); // ndof for x and y state
  } // end event loop

  // ============================================================
  // DRAW HISTOGRAMS
  // ============================================================

  // Canvas 1: Pull
  TCanvas *c_total = new TCanvas("c_total", "TOTAL", 1200, 1200);
  c_total->Divide(2, 3);
  c_total->cd(1);
  h_pull_y->Draw("");
  h_pull_y->Fit("gaus", "Q");
  c_total->cd(2);
  h_pull_x->Draw("");
  h_pull_x->Fit("gaus", "Q");

  // Canvas 2: Chi2
  c_total->cd(3);
  h_chi2_y->Draw("");
  c_total->cd(4);
  h_chi2_x->Draw("");

  // Canvas 3: Pvalue
  c_total->cd(5);
  h_pval->Draw();

  // Result Canvas
  TCanvas *c_result = new TCanvas("c_predict_momentum", "MOMENTUM", 1200, 500);
  c_result->Divide(2, 1);
  c_result->cd(1);
  h_momentum_results->Draw();
  c_result->cd(2);
  h_dp->Draw();
}

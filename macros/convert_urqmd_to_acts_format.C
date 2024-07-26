
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

void convert_urqmd_to_acts_format(){

  char buffer[300];
  float b,r0,rx,ry,rz,up0,upx,upy,upz,um;
  int evnr,ntracks,ityp,i2i3,ch,lcl,ncl,ior;

  TTree* t = new TTree("particles","particles");

  UInt_t event_id;
  vector<unsigned long> particle_id; 
  vector<int> particle_type; 
  vector<unsigned int> process;
  vector<float> vx;
  vector<float> vy;
  vector<float> vz;
  vector<float> vt;
  vector<float> px;
  vector<float> py;
  vector<float> pz;
  vector<float> m;
  vector<float> q;
  vector<float> eta;
  vector<float> phi;
  vector<float> pt;
  vector<float> p;
  vector<unsigned int>vertex_primary;
  vector<unsigned int>vertex_secondary;
  vector<unsigned int>particle;
  vector<unsigned int>generation;
  vector<unsigned int> sub_particle;
  t->Branch("particle_id",&particle_id);
  t->Branch("particle_type",&particle_type);
  t->Branch("process",&process);
  t->Branch("event_id",&event_id);
  t->Branch("vx",&vx);
  t->Branch("vy",&vy);
  t->Branch("vz",&vz);
  t->Branch("vt",&vt);
  t->Branch("px",&px);
  t->Branch("py",&py);  
  t->Branch("pz",&pz);
  t->Branch("m",&m);
  t->Branch("q",&q);
  t->Branch("eta",&eta);
  t->Branch("phi",&phi);  
  t->Branch("pt",&pt);
  t->Branch("p",&p);
  t->Branch("vertex_primary",&vertex_primary);
  t->Branch("vertex_secondary",&vertex_secondary);
  t->Branch("particle",&particle);
  t->Branch("generation",&generation);
  t->Branch("sub_particle",&sub_particle);
  

  FILE* ffile = fopen("urqmd.f14", "rb");
  TLorentzVector up;
  event_id = 0;
  while(!feof(ffile)) {
    particle_id.clear();
    particle_type.clear();
    process.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    vt.clear();
    px.clear();
    py.clear();
    pz.clear();
    m.clear();
    q.clear();
    eta.clear();
    phi.clear();
    pt.clear();
    p.clear();
    vertex_primary.clear();
    vertex_secondary.clear();
    particle.clear();
    generation.clear();
    sub_particle.clear();

    fgets(buffer, 300,ffile);  if (feof(ffile)) break;
    fgets(buffer, 300,ffile);
    fgets(buffer, 300,ffile);
    fgets(buffer, 36, ffile);
    fscanf(ffile, "%f", &b);
    fgets(buffer, 300, ffile);
    fgets(buffer, 300, ffile);
    fgets(buffer, 7, ffile);
    fscanf(ffile, "%d", &evnr);
    fgets(buffer, 300, ffile);
    for (int iline=0; iline<11; iline++)  { fgets(buffer, 300,ffile); }
    fscanf(ffile, "%d", &ntracks);
    fgets(buffer, 300, ffile);
    fgets(buffer, 300, ffile);
    int nCharged=1;
    for (int itrack = 0; itrack < ntracks; itrack++) {
      fgets(buffer, 300, ffile);
      double mm;
      sscanf(buffer,"%e %e %e %e %e %e %e %e %e %d %d %d  %d %d %d %*[^\n]",&r0,&rx,&ry,&rz,&up0,&upx,&upy,&upz,&um,&ityp,&i2i3,&ch,&lcl,&ncl,&ior);
      up.SetXYZM(upx,upy,upz,mm);
      float feta = up.Eta();
      //if (feta<1.5 || feta>3) continue;
      if      (ityp==   1 && ch== 1) { q.push_back( 1); particle_type.push_back(2212); mm = 0.93827208816; }
      else if (ityp== 101 && ch== 1) { q.push_back( 1); particle_type.push_back( 211); mm = 0.13957061;    }
      else if (ityp== 101 && ch==-1) { q.push_back(-1); particle_type.push_back(-211); mm = 0.13957061;    }
      else if (ityp== 106 && ch== 1) { q.push_back( 1); particle_type.push_back( 321); mm = 0.493677;      }
      else if (ityp==-106 && ch==-1) { q.push_back(-1); particle_type.push_back(-321); mm = 0.493677;      }
      else continue;
      
      uint64_t particleId = 0;
      particleId |= uint64_t(nCharged) << 52;
      particleId |= uint64_t(1) << 24;
      printf("%lu\n",particleId);
      particle_id.push_back(particleId);
      process.push_back(0);
      vx.push_back(0);
      vy.push_back(0);
      vz.push_back(0);
      vt.push_back(0);
      px.push_back(upx);
      py.push_back(upy);
      pz.push_back(upz);
      m.push_back(mm);
      eta.push_back(up.Eta());
      phi.push_back(up.Phi());
      pt.push_back(up.Pt());
      p.push_back(up.P());
      vertex_primary.push_back(nCharged);
      vertex_secondary.push_back(0);
      particle.push_back(1);
      generation.push_back(0);
      sub_particle.push_back(0);
      nCharged++;
    }
    t->Fill();
    event_id++;
  }
  fclose(ffile);

  TFile* fout = new TFile("urqmd.root","recreate");
  t->Write();
  fout->Close();
}

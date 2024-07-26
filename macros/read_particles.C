
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

void read_particles(){
  TFile* f = new TFile("urqmd.root");
  //TFile* f = new TFile("particles.root");
  f->ls();
  TTree* t = (TTree*) f->Get("particles");
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
  vector<unsigned long>* ptr_particle_id = &particle_id; 
  vector<int>* ptr_particle_type = &particle_type; 
  vector<unsigned int>* ptr_process = &process;
  vector<float>* ptr_vx = &vx;
  vector<float>* ptr_vy = &vy;
  vector<float>* ptr_vz = &vz;
  vector<float>* ptr_vt = &vt;
  vector<float>* ptr_px = &px;
  vector<float>* ptr_py = &py;
  vector<float>* ptr_pz = &pz;
  vector<float>* ptr_m = &m;
  vector<float>* ptr_q = &q;
  vector<float>* ptr_eta = &eta;
  vector<float>* ptr_phi = &phi;
  vector<float>* ptr_pt = &pt;
  vector<float>* ptr_p = &p;
  vector<unsigned int>* ptr_vertex_primary = &vertex_primary;
  vector<unsigned int>* ptr_vertex_secondary = &vertex_secondary;
  vector<unsigned int>* ptr_particle = &particle;
  vector<unsigned int>* ptr_generation = &generation;
  vector<unsigned int>* ptr_sub_particle = &sub_particle;

  t->SetBranchAddress("event_id",&event_id);
  t->SetBranchAddress("particle_id",&ptr_particle_id);
  t->SetBranchAddress("particle_type",&ptr_particle_type);
  t->SetBranchAddress("process",&ptr_process);
  t->SetBranchAddress("vx",&ptr_vx);
  t->SetBranchAddress("vy",&ptr_vy);
  t->SetBranchAddress("vz",&ptr_vz);
  t->SetBranchAddress("vt",&ptr_vt);
  t->SetBranchAddress("px",&ptr_px);
  t->SetBranchAddress("py",&ptr_py);  
  t->SetBranchAddress("pz",&ptr_pz);
  t->SetBranchAddress("m",&ptr_m);
  t->SetBranchAddress("q",&ptr_q);
  t->SetBranchAddress("eta",&ptr_eta);
  t->SetBranchAddress("phi",&ptr_phi);  
  t->SetBranchAddress("pt",&ptr_pt);
  t->SetBranchAddress("p",&ptr_p);
  t->SetBranchAddress("vertex_primary",&ptr_vertex_primary);
  t->SetBranchAddress("vertex_secondary",&ptr_vertex_secondary);
  t->SetBranchAddress("particle",&ptr_particle);
  t->SetBranchAddress("generation",&ptr_generation);
  t->SetBranchAddress("sub_particle",&ptr_sub_particle);
  for (uint64_t ev=0;ev<t->GetEntries();ev++){
    t->GetEntry(ev);
    for (int i=0;i<particle_id.size();i++){
      printf("%lu\n",particle_id[i]);
      printf("%u\n",particle_type[i]);
      printf("%u\n",process[i]);
      printf("%u\n",particle[i]);
      printf("%u\n",generation[i]);
      printf("%u\n",sub_particle[i]);
    }
  }
}

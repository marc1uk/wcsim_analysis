{
#include <deque>
#include <iostream>
#include <sstream>
#include <string>
// for reading csv files 
typedef std::deque <std::string> record_t;
typedef std::deque <record_t>    table_t;
std::istream& operator >> ( std::istream& ins, table_t& table ){
  std::string s;
  table.clear();
  while (std::getline( ins, s ))
    {
    std::istringstream ss( s );
    record_t           record;
    std::string        field;
    bool               final = true;
    while (std::getline( ss, field, ';' ))
      {
      record.push_back( field );
      final = ss.eof();
      }
    if (!final)
      record.push_back( std::string() );
    table.push_back( record );
    }
  return ins;
}

TFile* f = new TFile("HVlogs.root","RECREATE");
TTree* t= new TTree("hvlogs","Log of HV readings");

const char* inputdir="/media/marc/SAMSUNG/HV logs/Logs/Channel_Log_401_08_00/";
const char* ext = ".csv";
// "Monitoring_Log_1-5-8.csv";
TString nextfilepattern;
TSystemDirectory dir(inputdir, inputdir);
TList *subfolders = dir.GetListOfFiles(); // despite name returns files and folders
bool firstfile=true; // first file must also specify the branch structure
//const char* branchstructure ="On/b:Status/C:Chan_Name/C:Crate/b:Card_Type/b:Card/b:Channel/s:Vswmax/i:Vmon/i:Vset/i:Imon/i:Iset/i:Timestamp/C:Flagged/b:
// sample file line:
// 1,,TNK_0-6-0-0,1,,5,8,0.00,1700.80,1700.00,62.00uA,120.00uA,21/01/2017 03:17:33,0
// written from LabVIEW as an array of:
// ON/OFF | Status String | Channel Name | Crate | Card Type (enum) | Card | Channel | Vswmax |
// Vmon | Vset | Imon | Iset | Timestamp String | Channel Flagged
// Voltages's are in V, Currents have units specified using SI suffix
//C:a character string terminated by the 0 character
//B:an 8 bit signed integer
//b:an 8 bit unsigned integer
//S:a 16 bit signed integer
//s:a 16 bit unsigned integer
//I:a 32 bit signed integer
//i:a 32 bit unsigned integer
//L:a 64 bit signed integer
//l:a 64 bit unsigned integer
//F:a 32 bit floating point
//D:a 64 bit floating point
//It is also possible to manually specify the number of bytes to use instead of using the above: e.g. "ntrack/I2" will create a 2-byte (16-bit) int rather than a 32-bit one.

bool poweron, flagged;
std::string status, channel_name, timestampstring;
int crate, channel_polarity, card, channel;
double Vswmax, Vmon, Vset, Imon, Iset;
TTimestamp timestamp_ts;

// ON/OFF | Status String | Channel Name | Crate | Card Type (enum) | Card | Channel | Vswmax |
// Vmon | Vset | Imon | Iset | Timestamp String | Channel Flagged
t->Branch("Power",&poweron);
t->Branch("Status",&status);
t->Branch("Channel_Name",&channel_name);
t->Branch("Crate",&crate);
t->Branch("Polarity",&channel_polarity);
t->Branch("Card",&card);
t->Branch("Channel",&channel);
t->Branch("Vswmax",&Vswmax);
t->Branch("Vmon",&Vmon);
t->Branch("Vset",&Vset);
t->Branch("Imon",&Imon);
t->Branch("Iset",&Iset);
t->Branch("Timestamp_String");
t->Branch("Timestamp",&timestamp_ts);
t->Branch("Flagged",&flagged);

if(subfolders&&(subfolders->GetEntries()>2)) { // always at least 2 subfolders: '.' and '..'
  cout<<"looping over subfolders"<<endl;
  TSystemDirectory *subfolder;
  TIter nextsf(subfolders);
  while ((subfolder=(TSystemDirectory*)nextsf())) {
    if (subfolder->IsDirectory()){
      TString sfname = subfolder->GetName();
      if(sfname=="."||sfname==".."){ continue; }
      cout<<"found subfolder '"<<sfname<<"'"<<endl;
      /*
      // don't actually need this information
      TList *files = subfolder->GetListOfFiles();
      TSystemFile *file;
      TIter nextf(files);
      while (file=(TSystemFile*)nextf()) {
        TString fname;
        fname = file->GetName();
        if(fname.EndsWith(ext)) { int i=0; }
      }
      */
      nextfilepattern = inputdir + sfname + "/Monitoring_Log_*"; // XXX
      cout<<"adding "<<nextfilepattern<<endl;
      //t->Add(nextfilepattern.Data());
      //t->ReadFile(nextfilepattern.Data());
      FillTreeWithFile(nextfilepattern,t);
      // cf. t->Add("/home/marc/anniegpvm/stats10k/stats_1_79283/wcsim_*");
    } // end if directory
  } // end loop over subdirectories
} // end if subfolders exist.
if(1 /*t->GetEntriesFast()==0*/) {  // add files in parent directory
  // can't '+' const chars, convert to string
  nextfilepattern = std::string(inputdir) + std::string("/wcsim_0.*.root"); // XXX /wcsim_0.*
  cout<<"adding "<<nextfilepattern<<endl;
  t->Add(nextfilepattern.Data());
}
cout<<"Loaded "<<t->GetEntriesFast()<<" entries"<<endl;  // fast doesn't work with a tchain

void FillTreeWithFile(std::string filepath, TTree* t){
  get_file_contents(filepath.c_str());
}
bool poweron, flagged;
std::string status, channel_name, timestampstring;
int crate, channel_polarity, card, channel;
double Vswmax, Vmon, Vset, Imon, Iset;
TTimestamp timestamp_ts;

#include <fstream>
#include <string>
#include <cerrno>

std::string get_file_contents(const char *filename){
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in){
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
  throw(errno);
}



/*
Error file format:
  Timestamp String | Erro Code | Error Source | Screencap String | Serial Read String | Cursor Position

Event log format:
Timestamp | Event Type | Info

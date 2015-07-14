//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVectorD.h>

#include "Messenger/Messenger.h"
#include "Utils/XmlParserUtils.h"

using std::ostringstream;

//____________________________________________________________________________
xmlNodePtr genie::utils::xml::FindNode(xmlDocPtr xml_doc, string node_path)
{
  xmlNodePtr root_node = xmlDocGetRootElement(xml_doc);
  if(root_node==NULL) {
     LOG("XML", pERROR) << "Null root XML node";
     return NULL;
  }

  vector<string> node_path_vec = genie::utils::str::Split(node_path,"/");

  unsigned int ndepth = node_path_vec.size();
  unsigned int idepth = 0;

  xmlNodePtr curr_node = root_node; 

  while (curr_node != NULL) {
    if( (!xmlStrcmp(
       curr_node->name, (const xmlChar *) node_path_vec[idepth].c_str())) ) 
    {
        idepth++;
        if(idepth == ndepth) {
          return curr_node;
        } else {
          curr_node = curr_node->xmlChildrenNode;
        }
    }
    curr_node = curr_node->next;
  }

  xmlFree(curr_node);
  return NULL;
}
//____________________________________________________________________________
bool genie::utils::xml::GetBool(xmlDocPtr xml_doc, string node_path)
{
  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return false;
  }
  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );

  if(content == "true" || 
     content == "TRUE" || 
     content == "True" || 
     content == "on"   ||
     content == "ON"   ||
     content == "On"   ||
     content == "1"    || 
     content == "I") return true;

  if(content == "false" || 
     content == "FALSE" || 
     content == "False" ||
     content == "off"   ||
     content == "OFF"   ||
     content == "Off"   ||
     content == "0"     || 
     content == "O") return false;

  LOG("XML", pERROR) 
       << "Could not interpret content (" << content 
       << ") found at node_path: " << node_path << " as a boolean!";
  return false;
}
//____________________________________________________________________________
int genie::utils::xml::GetInt(xmlDocPtr xml_doc, string node_path)
{
  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return -999999;
  }
  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );
  int value = atoi(content.c_str());
  return value;
}
//____________________________________________________________________________
vector<int> genie::utils::xml::GetIntArray(xmlDocPtr xml_doc, string node_path)
{
  vector<int> array;

  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return array;
  }

  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );
  
  vector<string> str_tokens = genie::utils::str::Split(content, ",");
  vector<string>::const_iterator str_tokens_iter = str_tokens.begin();
  for( ; str_tokens_iter != str_tokens.end(); ++str_tokens_iter) {
     string token = genie::utils::str::TrimSpaces(*str_tokens_iter);
     if (not token.size() ) continue;
     array.push_back( atoi(token.c_str()) );
  }
  return array;
}
//____________________________________________________________________________
double genie::utils::xml::GetDouble(xmlDocPtr xml_doc, string node_path)
{
  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return -999999;
  }
  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );
  double value = atof(content.c_str());
  return value;
}
//____________________________________________________________________________
vector<double> 
  genie::utils::xml::GetDoubleArray(xmlDocPtr xml_doc, string node_path)
{
  vector<double> array;

  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return array;
  }

  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );
  
  vector<string> str_tokens = genie::utils::str::Split(content, ",");
  vector<string>::const_iterator str_tokens_iter = str_tokens.begin();
  for( ; str_tokens_iter != str_tokens.end(); ++str_tokens_iter) {
     string token = genie::utils::str::TrimSpaces(*str_tokens_iter);
     if (not token.size() ) continue;
     array.push_back( atof(token.c_str()) );
  }
  return array;
}
//____________________________________________________________________________
string genie::utils::xml::GetString(xmlDocPtr xml_doc, string node_path)
{
  xmlNodePtr node = genie::utils::xml::FindNode(xml_doc, node_path);
  if(node==NULL) {
    return "";
  }
  string content = genie::utils::xml::TrimSpaces( 
      xmlNodeListGetString(xml_doc, node->xmlChildrenNode, 1) );
  return content;
}
//____________________________________________________________________________
string genie::utils::xml::GetROOTFileName(xmlDocPtr xml_doc, string node_path)
{
  return genie::utils::xml::GetString(xml_doc, node_path+"/filename");
}
//____________________________________________________________________________
string genie::utils::xml::GetROOTObjName (xmlDocPtr xml_doc, string node_path)
{
  return genie::utils::xml::GetString(xml_doc, node_path+"/objname");
}
//____________________________________________________________________________
string genie::utils::xml::GetROOTObjType (xmlDocPtr xml_doc, string node_path)
{
  return genie::utils::xml::GetString(xml_doc, node_path+"/objtype");
}
//____________________________________________________________________________
TFile * genie::utils::xml::GetTFile(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  LOG("XML", pINFO) << "Reading info at XML node node_path: " << node_path;

  string filename = genie::utils::xml::GetROOTFileName(xml_doc,node_path);

  string used_base_dir = base_dir;
  if(base_dir.find("<env>") != string::npos) {
    used_base_dir = gSystem->Getenv("GENIE");
  }

  const char * full_filename = 
     Form("%s/%s", used_base_dir.c_str(), filename.c_str());
  TFile * file = new TFile(full_filename,"read");
  if(!file) {
    LOG("XML",pERROR) << "No file: " << full_filename << " found";
    return 0;
  }
  if(file->IsZombie()) {
    LOG("XML",pERROR) << "File is a zombie: " << full_filename;
    return 0;
  }

  return file;
}
//____________________________________________________________________________
TH1F * genie::utils::xml::GetTH1F(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  TFile * file = genie::utils::xml::GetTFile(xml_doc,node_path,base_dir);
  if(!file) return 0;

  string objname = genie::utils::xml::GetROOTObjName(xml_doc,node_path);
  TH1F * h = (TH1F*)file->Get(objname.c_str());
  if(!h) {
    LOG("XML",pERROR) << "No " << objname;
    file->Close();
    delete file;
    return 0;
  }

  TH1F * hcopy = new TH1F(*h);
  hcopy->SetDirectory(0);
  file->Close();
  delete file;

  return hcopy;
}
//____________________________________________________________________________
TH1D * genie::utils::xml::GetTH1D(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  TFile * file = genie::utils::xml::GetTFile(xml_doc,node_path,base_dir);
  if(!file) return 0;

  string objname = genie::utils::xml::GetROOTObjName(xml_doc,node_path);
  TH1D * h = (TH1D*)file->Get(objname.c_str());
  if(!h) {
    LOG("XML",pERROR) << "No " << objname;
    file->Close();
    delete file;
    return 0;
  }

  TH1D * hcopy = new TH1D(*h);
  hcopy->SetDirectory(0);
  file->Close();
  delete file;

  return hcopy;
}
//____________________________________________________________________________
TH2D * genie::utils::xml::GetTH2D(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  TFile * file = genie::utils::xml::GetTFile(xml_doc,node_path,base_dir);
  if(!file) return 0;

  string objname = genie::utils::xml::GetROOTObjName(xml_doc,node_path);
  TH2D * h = (TH2D*)file->Get(objname.c_str());
  if(!h) {
    LOG("XML",pERROR) << "No " << objname;
    file->Close();
    delete file;
    return 0;
  }

  TH2D * hcopy = new TH2D(*h);
  hcopy->SetDirectory(0);
  file->Close();
  delete file;

  return hcopy;
}
//____________________________________________________________________________
TVectorD * genie::utils::xml::GetTVectorD(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  TFile * file = genie::utils::xml::GetTFile(xml_doc,node_path,base_dir);
  if(!file) return 0;

  string objname = genie::utils::xml::GetROOTObjName(xml_doc,node_path);
  TVectorD * vec = (TVectorD*)file->Get(objname.c_str());
  if(!vec) {
    LOG("XML",pERROR) << "No " << objname;
    file->Close();
    delete file;
    return 0;
  }

  TVectorD * veccopy = new TVectorD(*vec);
  file->Close();
  delete file;

  return veccopy;
}
//____________________________________________________________________________
/*
TMatrixDSym * genie::utils::xml::GetTMatrixDSym(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  return 0;
}
//____________________________________________________________________________
TMatrixD * genie::utils::xml::GetTMatrixD(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  return 0;
}
//____________________________________________________________________________
TSpline3 * genie::utils::xml::GetTSpline3(
  xmlDocPtr xml_doc, string node_path, string base_dir)
{
  return 0;
}
//____________________________________________________________________________
*/

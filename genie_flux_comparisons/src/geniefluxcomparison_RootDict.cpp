// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIgeniefluxcomparison_RootDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "ColourWheel.hh"
#include "MRDspecs.hh"

// Header files passed via #pragma extra_include

namespace MRDSpecs {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *MRDSpecs_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("MRDSpecs", 0 /*version*/, "MRDspecs.hh", 12,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &MRDSpecs_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *MRDSpecs_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *ColourWheel_Dictionary();
   static void ColourWheel_TClassManip(TClass*);
   static void *new_ColourWheel(void *p = 0);
   static void *newArray_ColourWheel(Long_t size, void *p);
   static void delete_ColourWheel(void *p);
   static void deleteArray_ColourWheel(void *p);
   static void destruct_ColourWheel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ColourWheel*)
   {
      ::ColourWheel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ColourWheel));
      static ::ROOT::TGenericClassInfo 
         instance("ColourWheel", "ColourWheel.hh", 3,
                  typeid(::ColourWheel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ColourWheel_Dictionary, isa_proxy, 4,
                  sizeof(::ColourWheel) );
      instance.SetNew(&new_ColourWheel);
      instance.SetNewArray(&newArray_ColourWheel);
      instance.SetDelete(&delete_ColourWheel);
      instance.SetDeleteArray(&deleteArray_ColourWheel);
      instance.SetDestructor(&destruct_ColourWheel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ColourWheel*)
   {
      return GenerateInitInstanceLocal((::ColourWheel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ColourWheel*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ColourWheel_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ColourWheel*)0x0)->GetClass();
      ColourWheel_TClassManip(theClass);
   return theClass;
   }

   static void ColourWheel_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","include/ColourWheel.hh");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_ColourWheel(void *p) {
      return  p ? new(p) ::ColourWheel : new ::ColourWheel;
   }
   static void *newArray_ColourWheel(Long_t nElements, void *p) {
      return p ? new(p) ::ColourWheel[nElements] : new ::ColourWheel[nElements];
   }
   // Wrapper around operator delete
   static void delete_ColourWheel(void *p) {
      delete ((::ColourWheel*)p);
   }
   static void deleteArray_ColourWheel(void *p) {
      delete [] ((::ColourWheel*)p);
   }
   static void destruct_ColourWheel(void *p) {
      typedef ::ColourWheel current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ColourWheel

namespace {
  void TriggerDictionaryInitialization_geniefluxcomparison_RootDict_Impl() {
    static const char* headers[] = {
"ColourWheel.hh",
"MRDspecs.hh",
0
    };
    static const char* includePaths[] = {
"/grid/fermiapp/products/larsoft/root/v6_06_08/Linux64bit+2.6-2.12-e10-nu-debug/include",
"/grid/fermiapp/products/larsoft/genie/v2_12_0a/Linux64bit+2.6-2.12-e10-r6-debug/GENIE_R2120/src",
"/grid/fermiapp/products/larsoft/root/v6_06_08/Linux64bit+2.6-2.12-e10-nu-debug/include",
"/annie/app/users/moflaher/wcsim/root_work/genie_flux_comparisons/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "geniefluxcomparison_RootDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(file_name@@@include/ColourWheel.hh)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$ColourWheel.hh")))  ColourWheel;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "geniefluxcomparison_RootDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "ColourWheel.hh"
#include "MRDspecs.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"ColourWheel", payloadCode, "@",
"ColourWheel::colournames", payloadCode, "@",
"ColourWheel::colours", payloadCode, "@",
"MRDSpecs::MRDPMTExposeHeight", payloadCode, "@",
"MRDSpecs::MRDPMTRadius", payloadCode, "@",
"MRDSpecs::MRD_depth", payloadCode, "@",
"MRDSpecs::MRD_end", payloadCode, "@",
"MRDSpecs::MRD_height", payloadCode, "@",
"MRDSpecs::MRD_layer2", payloadCode, "@",
"MRDSpecs::MRD_start", payloadCode, "@",
"MRDSpecs::MRD_steel_height", payloadCode, "@",
"MRDSpecs::MRD_steel_width", payloadCode, "@",
"MRDSpecs::MRD_width", payloadCode, "@",
"MRDSpecs::alufullxlen", payloadCode, "@",
"MRDSpecs::alufullxthickness", payloadCode, "@",
"MRDSpecs::alufullylen", payloadCode, "@",
"MRDSpecs::alufullythickness", payloadCode, "@",
"MRDSpecs::alufullzlen", payloadCode, "@",
"MRDSpecs::alusteelgap", payloadCode, "@",
"MRDSpecs::fidcutradius", payloadCode, "@",
"MRDSpecs::fidcuty", payloadCode, "@",
"MRDSpecs::fidcutz", payloadCode, "@",
"MRDSpecs::heights", payloadCode, "@",
"MRDSpecs::layergap", payloadCode, "@",
"MRDSpecs::layeroffsets", payloadCode, "@",
"MRDSpecs::maxheight", payloadCode, "@",
"MRDSpecs::maxwidth", payloadCode, "@",
"MRDSpecs::mrdZlen", payloadCode, "@",
"MRDSpecs::mrdpmtfullheight", payloadCode, "@",
"MRDSpecs::mrdscintlayers", payloadCode, "@",
"MRDSpecs::nothickness", payloadCode, "@",
"MRDSpecs::numalustructs", payloadCode, "@",
"MRDSpecs::numhpanels", payloadCode, "@",
"MRDSpecs::nummrdpmts", payloadCode, "@",
"MRDSpecs::numpaddlesperpanelh", payloadCode, "@",
"MRDSpecs::numpaddlesperpanelv", payloadCode, "@",
"MRDSpecs::numpanels", payloadCode, "@",
"MRDSpecs::numplates", payloadCode, "@",
"MRDSpecs::numvetopaddles", payloadCode, "@",
"MRDSpecs::numvpanels", payloadCode, "@",
"MRDSpecs::paddle_extentsx", payloadCode, "@",
"MRDSpecs::paddle_extentsy", payloadCode, "@",
"MRDSpecs::paddle_extentsz", payloadCode, "@",
"MRDSpecs::paddle_layers", payloadCode, "@",
"MRDSpecs::paddle_orientations", payloadCode, "@",
"MRDSpecs::paddle_originx", payloadCode, "@",
"MRDSpecs::paddle_originy", payloadCode, "@",
"MRDSpecs::paddle_originz", payloadCode, "@",
"MRDSpecs::scintalugap", payloadCode, "@",
"MRDSpecs::scintbordergap", payloadCode, "@",
"MRDSpecs::scintfullxlen", payloadCode, "@",
"MRDSpecs::scintfullzlen", payloadCode, "@",
"MRDSpecs::scinthfullylen", payloadCode, "@",
"MRDSpecs::scintlgfullheight", payloadCode, "@",
"MRDSpecs::scintlgfullwidth", payloadCode, "@",
"MRDSpecs::scinttapfullheight", payloadCode, "@",
"MRDSpecs::scinttapfullwidth", payloadCode, "@",
"MRDSpecs::scintvfullylen", payloadCode, "@",
"MRDSpecs::steelfullxlen", payloadCode, "@",
"MRDSpecs::steelfullylen", payloadCode, "@",
"MRDSpecs::steelfullzlen", payloadCode, "@",
"MRDSpecs::steelscintgap", payloadCode, "@",
"MRDSpecs::tank_halfheight", payloadCode, "@",
"MRDSpecs::tank_radius", payloadCode, "@",
"MRDSpecs::tank_start", payloadCode, "@",
"MRDSpecs::tank_yoffset", payloadCode, "@",
"MRDSpecs::tankouterRadius", payloadCode, "@",
"MRDSpecs::vetopaddlesperpanel", payloadCode, "@",
"MRDSpecs::widths", payloadCode, "@",
"MRDSpecs::windowheight", payloadCode, "@",
"MRDSpecs::windowwidth", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("geniefluxcomparison_RootDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_geniefluxcomparison_RootDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_geniefluxcomparison_RootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_geniefluxcomparison_RootDict() {
  TriggerDictionaryInitialization_geniefluxcomparison_RootDict_Impl();
}

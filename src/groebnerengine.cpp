#include "groebnerengine.h"
#include "printer.h"
#include "log.h"

GroebnerEngine *GroebnerEngine::list;



static GroebnerEngine *gfan,*singular,*default_;
static bool initialized;


GroebnerEngine::GroebnerEngine()
{
   next=list;
   list=this;
}


GroebnerEngine *GroebnerEngine::find(const char *name)
{
   GroebnerEngine *l=list;
   while(l)
      {
         if(std::string(l->name())==std::string(name))break;
         l=l->next;
      }
   return l;
}


static inline void GE_Init()
{
  if(!initialized)
    {
      GE_SetEngine("");
    }
}


bool GE_SetEngine(const char *name)
{
  singular=GroebnerEngine::find("singular");
  gfan=GroebnerEngine::find("gfan");
  GroebnerEngine *selected=GroebnerEngine::find(name);
  default_=gfan;
  if(singular)default_=singular;
  if(selected)default_=selected;
  initialized=true;
  assert(default_);
  log1 fprintf(stderr,"Groebner basis Engine being used: \"%s\".\n",default_->name());

  return selected;
}


PolynomialSet GE_groebnerBasis(PolynomialSet const &idealGenerators, TermOrder const &termOrder, bool autoreduce, bool allowSaturation)
{
  GE_Init();
  PolynomialSet ret(idealGenerators.getRing());
  bool success=true;

  //Should be put into a list instead

  if(default_)
    {
      ret=default_->groebnerBasis(success,idealGenerators,termOrder,autoreduce,allowSaturation);
      if(success)
	return ret;
    }

  if(singular)
    {
      ret=singular->groebnerBasis(success,idealGenerators,termOrder,autoreduce,allowSaturation);
      if(success)
	return ret;
    }
  if(gfan)
    {
      ret=gfan->groebnerBasis(success,idealGenerators,termOrder,autoreduce,allowSaturation);
      if(success)
	return ret;
    }
  assert(0);
  return ret;
}


PolynomialSet GE_autoReduce(PolynomialSet const &idealGenerators, class WeightRecipe const *recipe)
{
  GE_Init();
  PolynomialSet ret(idealGenerators.getRing());
  bool success=true;

  if(default_)
    {
      ret=default_->autoReduce(success,idealGenerators,recipe);
      if(success)
	return ret;
    }
  if(singular)
    {
      ret=singular->autoReduce(success,idealGenerators,recipe);
      if(success)
	return ret;
    }
  if(gfan)
    {
      ret=gfan->autoReduce(success,idealGenerators,recipe);
      if(success)
	return ret;
    }
  assert(0);
  return ret;
}

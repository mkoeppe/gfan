#pragma once
/*
 * gfanlib_memoryresource.h
 *
 *  Created on: 6 Jul 2018
 *      Author: anders
 */

#include <experimental/memory_resource>
#include <iostream>
#include <experimental/vector>
#include <assert.h>


/*In compiler explorer wih gcc8.1 -std=c++17 -O3
 *
 * #include <vector>

#include <experimental/memory_resource>
#include <experimental/vector>

class MyInt{
public:
	int a;
//	int *p=&a;
	//MyInt(){seta(a);}
	MyInt(int a_):a(a_){}
};


// Type your code here, or load an example.
void test() {
    				 MyInt t(0);
	//				 t.a=0;
				 for(int i=0;i<100000;i++)
					 {
						 std::experimental::pmr::vector<MyInt> A(100000,t/*	,&stackResource/);
						 //A.resize(100000,MyInt{});
					 }
			 }
 */

template<typename a>
//using std::experimental::pmr pmr;
using pmrvector=std::experimental::pmr::vector<a>;
typedef std::experimental::pmr::memory_resource MR; //Maybe this should be polymorphic_allocator<byte> instead - at least for the constructors using these
using std::experimental::pmr::get_default_resource;

//DELETE
class ResourceWrapper: public MR{
	MR *parent;
	int numberOfCurrentAllocations;
public:
	ResourceWrapper(MR *mr_):
		parent(mr_),
		numberOfCurrentAllocations(0)
	{
	}
	~ResourceWrapper()
	{
		assert(numberOfCurrentAllocations==0);
	}
	virtual void* do_allocate(size_t __bytes, size_t __alignment)
	{
		numberOfCurrentAllocations++;
		return parent->allocate(__bytes,__alignment);
	}
	virtual void do_deallocate(void* __p, size_t __bytes, size_t __alignment)
	{
		numberOfCurrentAllocations--;
		parent->deallocate(__p,__bytes,__alignment);
	}
	virtual bool do_is_equal(const memory_resource& __other) const noexcept
	{
		return this==&__other;
	}
};

class StackResource: public std::experimental::pmr::memory_resource{
  /*
   *
   * Layout of memory
   * The header is put after an allocation. It contains an offset from this header to the previous header.
   * Each header is 32 bit, with LSB being 1 if if the memory chunck before the header has been freed. The SuperHeader always has LSB 0;
   * Top header
   *
    Layout for memory:
    SuperHeader
    padding
    block
    padding
    Header
    padding
    block
    padding
    header  <---topHeader

	The header is 4 byte aligned, while the alignment of each block is specified at allocation.
   */
public:
  std::experimental::pmr::vector<char> mem;
  int topHeader;
  int nbytes;
  memory_resource* parentResource;
  memory_resource* fallBackResource;
  int64_t numberOfAllocatedBytes;
  int64_t numberOfAllocs;
  int depth;
  int maxDepth;
  int maxTopHeader;
  bool report;
  ~StackResource()
  {
	  if(report)
	  {
		  std::cerr<<"Number of allocs\t"<<numberOfAllocs<<"\n"
				  <<"MaxDepth\t\t"<<maxDepth<<"\n"
				  <<"MaxUsedBytes\t\t"<<maxTopHeader+4<<"\n";
	  }
  }
  int firstPosition(int firstPossible, int align)
  {
	  int minimumAlignment=4;
	  return (firstPossible+((align-1)|(minimumAlignment-1)))& ~((align-1)|(minimumAlignment-1));
  }
  StackResource(int nbytes_,
		memory_resource* parentResource_=std::experimental::pmr::new_delete_resource(),
		memory_resource* fallBackResource_=std::experimental::pmr::null_memory_resource()):
	topHeader(0),
    nbytes(nbytes_),
    mem(nbytes_,parentResource_),
	numberOfAllocatedBytes(0),
    parentResource(parentResource_),
    fallBackResource(fallBackResource_),
	numberOfAllocs(0),
	depth(0),
	maxDepth(0),
	maxTopHeader(0),
	report(0)//Enable reporting here. It would be great if this was an argument to the constructor.
  {
    assert(nbytes>3);
    *(int*)(&(mem[0]))=0;
  }
  int &intAt(int index)
  {
	  return (*(int*)&(mem[index]));
  }
  char stateOfIthLongWord(int i)
  {
	  int p=topHeader;
	  while(p)
	  {
		  int n=p+(intAt(p)&~1);
		  if(4*i==p)return 'H';
		  if(n<4*i && 4*i<p)return (intAt(p)&1?' ':'*');
		  p=n;
	  }
	  return 'I';
  }
  StackResource *dump()
  {
	  std::cerr<<"Resource dump:\n"<<"topHeader="<<topHeader<<"Bytes allocated"<<numberOfAllocatedBytes<<"\n";
	  for(int i=0;i<topHeader/4+1;i++)
	  {
		  std::cerr<<(int*)&(mem[i*4])<<":"<<i*4<<stateOfIthLongWord(i)<<":"<<"\t"<<((int*)&(mem[0]))[i]<<"\n";
	  }

	  {
		  int used=0;
		  int free=0;
		  int head=0;

		  for(int i=0;i<topHeader/4+1;i++)
		  {
			  switch(stateOfIthLongWord(i))
			  {
			  case 'I':head++;break;
			  case 'H':head++;break;
			  case '*':used++;break;
			  case ' ':free++;break;
			  }
		  }
		  int total=head+used+free;
		  std::cerr<<"Used:\t"<<int(used/float(total)*100)<<"%\n";
		  std::cerr<<"Head:\t"<<int(head/float(total)*100)<<"%\n";
		  std::cerr<<"Free:\t"<<int(free/float(total)*100)<<"%\n";
	  }
	  return this;
  }
  bool breakCondition(int address, int size)
  {
	  return false;
	  return (2876==address) && (size==168);
	  return (3048==address) && (size==168);
	  return (4252==address) && (size==936);
	  return (3992==address) && (size==256);
  }
  int adjustTopHeader()
  {
	  while(intAt(topHeader)&1)
	  {
		  topHeader+=(intAt(topHeader)&~1);
		  depth--;
	  }
	  return topHeader;
  }
  virtual void* do_allocate(size_t __bytes, size_t __alignment)
  {
	  numberOfAllocatedBytes+=__bytes;
//	  std::cerr<<"ALLOC\n";
	  numberOfAllocs++;
//  dump();
    if(__bytes>2000000000)
      return fallBackResource->allocate(__bytes,__alignment);
    adjustTopHeader();
    int ret=firstPosition(topHeader+4,__alignment);
    int nextfreeCandidate=firstPosition(ret+__bytes,4);
    if(nextfreeCandidate+4>nbytes)
    {
    	dump();assert(0);
      return fallBackResource->allocate(__bytes,__alignment);
    }
    intAt(nextfreeCandidate)=topHeader-nextfreeCandidate;
//    for(int i=ret;i<ret+__bytes;i++)mem[i]=0;//Clear
    topHeader=nextfreeCandidate;
//    std::cerr<<"ALLOC \t"<<__bytes<<" bytes \t"/*<<__alignment*/<<"AT\t"<<ret<<"\n";
//    dump();
    depth++;
    if(depth>maxDepth)maxDepth=depth;
    if(topHeader>maxTopHeader)maxTopHeader=topHeader;

    assert(!breakCondition(ret,__bytes));

    return &(mem[ret]);
  }

  virtual void do_deallocate(void* __p, size_t __bytes, size_t __alignment){
	  numberOfAllocatedBytes-=__bytes;
//    cerr<<"FREE\t"<<(char*)__p-(char*)(&(mem[0]))<<"\n";
//    dump();

    if(__p<(&(mem[0]))||__p>((&(mem[0]))+nbytes))
      fallBackResource->deallocate(__p,__bytes,__alignment);
    int nextHeader=firstPosition((((char*)__p)-((char*)&(mem[0])))+__bytes,4);
    intAt(nextHeader)|=1;
//    dump();
  }

  virtual bool
  do_is_equal(const memory_resource& __other) const noexcept
  {
    return this==&__other;
  }
};


std::string memoryResourceToString(const std::experimental::pmr::memory_resource *mr);

// DELETE?
class ResourceInvariant
{
	StackResource *sr;
	int bytesAllocatedAtStart;
	int topHeaderAtStart;
public:
	ResourceInvariant(MR *sr_):
		sr((StackResource*)sr_),
		bytesAllocatedAtStart((dynamic_cast<StackResource*>(sr_))->numberOfAllocatedBytes),
		topHeaderAtStart((dynamic_cast<StackResource*>(sr_))->adjustTopHeader())
	{
	}
	void assertIsInvariant()
	{
		if(bytesAllocatedAtStart!=sr->numberOfAllocatedBytes || topHeaderAtStart!=(dynamic_cast<StackResource*>(sr))->adjustTopHeader())
		{
			std::cerr<<"Number of allocated bytes at start:"<<bytesAllocatedAtStart<<"\n"
					<<"Number of allocated bytes at test:"<<sr->numberOfAllocatedBytes<<"\n"
					<<"Topheader at start:"<<topHeaderAtStart<<"\n"
					<<"Topheader at test:"<<(dynamic_cast<StackResource*>(sr))->adjustTopHeader()<<"\n";
			sr->dump();
		}
		assert(bytesAllocatedAtStart==sr->numberOfAllocatedBytes);
		assert(topHeaderAtStart==(dynamic_cast<StackResource*>(sr))->adjustTopHeader());
	}
	~ResourceInvariant()
	{
	}
};

extern StackResource stackResource;
extern StackResource stackResource2;

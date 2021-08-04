#ifndef PARALLELDEPEND_HPP_INCUDED
#define PARALLELDEPEND_HPP_INCUDED

#include <solve.hpp>
using namespace ngsolve;
#include "concurrentqueue.h"

typedef moodycamel::ConcurrentQueue<int> TQueue;
typedef moodycamel::ProducerToken TPToken;
typedef moodycamel::ConsumerToken TCToken;

static TQueue queue;

using namespace ngstd;

template <typename TFUNC>
void RunParallelDependency (const Table<int> & dag, TFUNC func)
{
  Array<atomic<int>> cnt_dep(dag.Size());

  for (auto & d : cnt_dep)
    d.store (0, memory_order_relaxed);


  ParallelFor (Range(dag),
	       [&] (int i)
	       {
		 for (int j : dag[i])
		   cnt_dep[j]++;
	       });


  Array<int> ready(dag.Size());
  ready.SetSize0();
  int num_final = 0;

  for (int j : Range(cnt_dep))
    {
      if (cnt_dep[j] == 0) ready.Append(j);
      if (dag[j].Size() == 0) num_final++;
    }


  /*
    while (ready.Size())
    {
    int size = ready.Size();
    int nr = ready[size-1];
    ready.SetSize(size-1);

    func(nr);

    for (int j : dag[nr])
    {
    cnt_dep[j]--;
    if (cnt_dep[j] == 0)
    ready.Append(j);
    }
    }
  */


  atomic<int> cnt_final(0);
  SharedLoop sl(Range(ready));

  Array< Vec<3> > timings(task_manager -> GetNumThreads());
  // double starttime = omp_get_wtime();


  task_manager -> CreateJob
    ([&] (const TaskInfo & ti)
     {
       TPToken ptoken(queue);
       TCToken ctoken(queue);

       for (int i : sl)
	 queue.enqueue (ptoken, ready[i]);

       while (1)
	 {
	   if (cnt_final >= num_final) break;

	   int nr;
	   if(!queue.try_dequeue_from_producer(ptoken, nr))
	     if(!queue.try_dequeue(ctoken, nr))
	       continue;

	   if (dag[nr].Size() == 0)
	     cnt_final++;

	   func(nr);

	   for (int j : dag[nr])
	     {
	       if (--cnt_dep[j] == 0)
		 queue.enqueue (ptoken, j);
	     }
	 }
     });

  /*
  // my own simple parall
  MyQueue<int> queue(dag.Size());

  for (int i : ready)
  queue.Push(i);

  task_manager -> CreateJob
  ([&] (const TaskInfo & ti)
  {
  while (1)
  {
  int nr;
  if (!queue.Pop(nr)) break;

  func(nr);

  for (int j : dag[nr])
  {
  if (--cnt_dep[j] == 0)
  queue.Push (j);
  }
  }
  ptqend[ti.task_nr] = omp_get_wtime();
  });
  */
}


#endif

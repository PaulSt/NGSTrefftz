TentPitchedSlab<DIM> tps
    = TentPitchedSlab<DIM> (ma); // collection of tents in timeslab

tps.PitchTents (adt, awavespeed); // adt = time slab height

Propagate (BaseVector &hu, LocalHeap &lh) // BaseVector ??? -> local storage
{
  AutoVector hu2 = hu.CreateVector ();

  RunParallelDependency (tps.tent_dependency, [&] (int i) {
    LocalHeap slh = lh.Split (); // split to threads
    Tent &tent = *tps.tents[i];
  });
}

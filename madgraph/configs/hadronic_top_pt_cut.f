C     Custom MadGraph cut for the current semileptonic ttbar process.
C     External particles 3, 4, and 5 are the hadronic top daughters (b j j).
C     MIN_HADRONIC_TOP_PT is read from the run card via run.inc.
      logical function dummy_cuts(P)
      implicit none
      include 'nexternal.inc'
      include 'maxparticles.inc'
      include 'run.inc'

      double precision P(0:3,nexternal)
      double precision px_had_top, py_had_top, pt_had_top

      dummy_cuts = .true.

      px_had_top = P(1,3) + P(1,4) + P(1,5)
      py_had_top = P(2,3) + P(2,4) + P(2,5)
      pt_had_top = dsqrt(px_had_top**2 + py_had_top**2)

      if (pt_had_top .lt. min_hadronic_top_pt) then
         dummy_cuts = .false.
         return
      endif

      return
      end

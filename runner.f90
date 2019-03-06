program kd_tester

   use kd_tree, only : kdnode
   use kd_tree, only : kd_insert, kd_construct, kd_search, kd_remove
   use kd_tree, only : kd_free
   use kd_tree, only : kd_find_min
   use getoptf

   implicit none

   integer :: argc
   character (len=255) :: argv
   character :: c

   integer :: npoints
   integer :: ndims

   integer :: lrange
   integer :: urange

   integer :: ntests
   integer :: npoints_l
   integer :: npoints_u
   integer :: ndims_l
   integer :: ndims_u

   namelist /testVars/ npoints, ndims, lrange, urange, ntests, npoints_l, npoints_u, ndims_l, ndims_u

   open(42,file='namelist.input', status='old',form='formatted',action='read')
   read(42, testVars)
   close(42)

   argc = command_argument_count()
   call get_command(argv)

   do while ( getopt(argc, argv, c, "12345d") )
      select case (c)
         case ('2')
            write(0,*) "Launching test 2"
            call test2(npoints, ndims, lrange, urange)
         case ('3')
            write(0,*) "Launching test 3"
            call test3(npoints, ndims, lrange, urange)
         case ('4')
            write(0,*) "Launching test 4"
            call test4(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('5')
            write(0,*) "Launching test 5"
            call test5(ntests, npoints, ndims, lrange, urange)
         case ('6')
            write(0,*) "Launching test 6"
            call test6(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('m')
            write(0,*) "Launching Minimum Test"
            call test_min(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('r')
            write(0,*) "Launching remove Test"
            call test_remove(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('d')
            write(0,*) "Launching remove Test"
            call test_distance2(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('i')
            write(0,*) "Launching insert test"
            call test_insert(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
      endselect
   enddo

contains

subroutine test_distance2(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none
   
   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(ndims) :: min_point 
   integer :: i, j, tests
   real :: r_ndims, r_npoints
   integer :: ndims, npoints
   real :: min_d, dis, bf_min
   real, dimension(:), allocatable :: minimum
   real, dimension(:), allocatable :: result, actuall_min_point

   do tests = 1, ntests, 1

      call random_number(r_ndims)
      call random_number(r_npoints)

      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) ""
      write(0,*) "Test number: ", tests
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))
      allocate(minimum(ndims))
      allocate(result(ndims))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      tree => kd_construct(arry1(:,:))

      do i = 1, size(arry2, dim=2)

         min_d = huge(min_d)
         call kd_search(tree, arry2(:,i), result, min_d)

         bf_min = huge(bf_min)
         dis = huge(dis)

         do j = 1, size(arry1, dim=2), 1
            dis = sum((arry2(:,i) - arry1(:,j))**2)
            if (dis < bf_min) then
               bf_min = dis
               actuall_min_point = arry1(:,j)
            endif
         end do

         if (bf_min /= min_d) then
            call print_tree(tree, 0, 0)
            write(0,*) "We did not find the nearest neighbor :("
            write(0,*) "Search Distance: ", min_d
            write(0,*) "Actual Distance: ", bf_min
            write(0,*) "We were searching for: "
            write(0,*) arry2(:,i)
            write(0,*) "KD Search found: "
            write(0,*) result(:)
            write(0,*) "BF actual nearest neighbor: "
            write(0,*) actuall_min_point(:)
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      deallocate(minimum)
      deallocate(result)
      call kd_free(tree)
  enddo

end subroutine test_distance2

subroutine test_insert(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
   implicit none
   
   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(:), pointer :: result
   integer :: i, j, tests
   real :: r_ndims, r_npoints
   integer :: ndims2, npoints2
   real :: min_d, dis, bf_min

   do tests = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) ""
      write(0,*) "Test number: ", tests
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))
      allocate(result(ndims))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      do i = 1, size(arry2, dim=2)
         call kd_insert(tree, arry1(:, i))

         min_d = huge(min_d)
         call kd_search(tree, arry1(:,i), result, min_d)
         if ( .NOT. all(arry1(:,i) == result(:))) then
            write(0,*) "We could not find the node we just inserted.. :("
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      deallocate(result)
      call kd_free(tree)
  enddo

end subroutine test_insert

subroutine test2(npoints, ndims, lrange, urange)

   implicit none
   integer, intent(in), value :: npoints, ndims
   integer, intent(in), value :: lrange, urange

   type(kdnode), pointer :: rTree => null()
   real, dimension(ndims, npoints) :: arry1, arry2
   real, dimension(ndims) :: result
   real :: min_d

   integer :: i, r, k, n

   min_d = huge(min_d)

   write(0,*) "Real test: "
   write(0,*) "Num Points: ", npoints
   write(0,*) "Num Dims: ", ndims

   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange
   
   rTree => kd_construct(arry1)

   write(0,*) ""

   if (associated(rTree)) then
      write(0,*) "Real Tree Node: ", rTree % point

      if (associated(rTree % left)) then
         write(0,*) "Left tree was associated"
         write(0,*) rTree % left % point
      else
         write(0,*) "Left tree was NOT associated"
      endif

      if (associated(rTree % right)) then
         write(0,*) "Right tree was associated"
         write(0,*) rTree % right % point
      else
         write(0,*) "Right tree was NOT associated"
      endif
   else
      write(0,*) "rTree was not associated :'("
   endif

   write(0,*) ""
   write(0,*) "Searching Test!"
   write(0,*) ""
   write(0,*) "Searching all points currently in the tree"
   write(0,*) ""

   do i = 1, size(arry1, dim=2)
      min_d = huge(min_d)
      call kd_search(rTree, arry1(:, i), result, min_d)
      if (.NOT. (all(result(:) == arry1(:,i)) .AND. min_d == 0.0)) then
         write(0,*) "That point wasn't in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      endif
   enddo

   write(0,*) ""
   write(0,*) "Searching for randomly created points"
   write(0,*) ""

   do i = 1, size(arry2, dim=2)
      min_d = huge(min_d)
      call kd_search(rTree, arry2(:, i), result, min_d)
      if (.NOT. (all(result(:) == arry2(:,i)) .AND. min_d == 0.0)) then
         write(0,*) "That point wasn't in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      endif
   enddo

   call kd_free(rTree)

end subroutine test2

subroutine test3(npoints, ndims, lrange, urange)

   use iso_c_binding, only : c_int

   implicit none

   interface
       subroutine timer_start(timer_id) bind(C)
          use iso_c_binding, only : c_int
          integer (c_int), intent(in), value :: timer_id
       end subroutine timer_start

       subroutine timer_stop(timer_id, sec, nsec) bind(C)
          use iso_c_binding, only : c_int
          integer (c_int), intent(in), value :: timer_id
          integer (c_int), intent(out) :: sec, nsec
       end subroutine timer_stop
   end interface
 
   integer (c_int) :: timer_id, sec, nsec 
   integer, intent(in), value :: npoints, ndims
   integer, intent(in), value :: lrange, urange

   type(kdnode), pointer :: rTree => null()
   real, dimension(ndims, npoints) :: arry1, arry2
   real, dimension(ndims) :: result
   real :: min_d

   integer :: i, r, k, n

   min_d = huge(min_d)

   write(0,*) "Timing test. Size: ", npoints, " dims: ", ndims
   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

   write(0,*) "Constructing tree..." 
   call timer_start(0)
   rTree => kd_construct(arry1) 
   call timer_stop(0, sec, nsec)
   write(0,*) "Tree constructed... time: ", sec, nsec

   
   write(0,*) "Searching for points that are within the tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call kd_search(rTree, arry1(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Searching for points that are NOT within the tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call kd_search(rTree, arry2(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Freeing rTree"
   call timer_start(0)
   call kd_free(rTree)
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "=============================================="
   write(0,*) "Testing timing for inserting nodes"
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call kd_insert(rTree, arry1(:, i))
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "Tree insertion creation time: ", sec, nsec

   write(0,*) "Searching for points that are within the INSERTED tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call kd_search(rTree, arry1(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Searching for points that are NOT within the INSERTED tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call kd_search(rTree, arry2(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Freeing rTree"
   call timer_start(0)
   call kd_free(rTree)
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   
end subroutine test3

function search_tree(tree, point)

   implicit none

   type(kdnode), pointer, intent(in) :: tree
   real, dimension(:), intent(in) :: point
   real, dimension(size(point, dim=1)) :: result

   logical search_tree
   real :: min_d

   search_tree = .TRUE.

   min_d = huge(min_d)
   !write(0,*) "New Search", point
   call kd_search(tree, point, result, min_d)
   if (.NOT. (all(result(:) == point(:)) .AND. min_d == 0.0)) then
      write(0,*) point(:), " was not in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      search_tree = .FALSE.
      return
   endif

end function search_tree

subroutine test4(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   real, dimension(:,:), pointer :: arry1, arry2
   type(kdnode), pointer :: tree => null()

   real :: r_ndims, r_npoints
   integer :: ndims, npoints

   integer :: i, j

   write(0,*) "Running Random tests"
   write(0,*) " running ", ntests, " tests"
   write(0,*) ""

   do i = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", i
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      tree => kd_construct(arry1)
      
      !write(0,*) "Arry1: ", arry1

      do j = 1, size(arry1, dim=2), 1
         if( .NOT. search_tree(tree, arry1(:,j))) then
            write(0,*) "This point was not found"
            call print_tree(tree, 0, 0)
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      call kd_free(tree)

   enddo

end subroutine test4

subroutine test5(ntests, npoints, ndims, lrange, urange)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: npoints
   integer, intent(in), value :: ndims
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   integer :: i, j
   real :: r_ndims, r_npoints
   integer :: ndims2, npoints2

   allocate(arry1(ndims, npoints))
   allocate(arry2(ndims, npoints))

   call random_number(arry1(:,:))
   call random_number(arry2(:,:))

   arry1(:,:) = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2(:,:) = (arry2(:,:) * (urange + 1 - lrange)) + lrange

   write(0,*) "Number of points: ", npoints, " Number of dims: ", ndims

   do i = 1, size(arry1, dim=2), 1
      call kd_insert(tree, arry1(:,i))
   enddo

   write(0,*) "All ", npoints, "were added into the tree!"

   do i = 1, size(arry1, dim=2), 1
      if( .NOT. search_tree(tree, arry1(:,i))) then
         write(0,*) "This point was not found! But it should have been!!"
         stop
      endif
   enddo

   write(0,*) "All points that were inserted to make the tree were succesfully found within the tree!!"
   call kd_free(tree)

end subroutine test5

subroutine test6(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   real, dimension(:,:), pointer :: arry1, arry2
   type(kdnode), pointer :: tree => null()

   real :: r_ndims, r_npoints
   integer :: ndims, npoints

   integer :: i, j

   write(0,*) "Running Random tests with insert"
   write(0,*) " running ", ntests, " tests"
   write(0,*) ""

   do i = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", i
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      do j = 1, npoints, 1
         call kd_insert(tree, arry1(:, j))
      enddo
      
      !write(0,*) "Arry1: ", arry1

      do j = 1, size(arry1, dim=2), 1
         if( .NOT. search_tree(tree, arry1(:,j))) then
            write(0,*) "This point was not found"
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      call kd_free(tree)

   enddo

end subroutine test6

subroutine test_min(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(ndims) :: min_point 
   real, dimension(2,7) :: arry3
   real :: r_ndims, r_npoints, r_test_dim
   integer :: test_dim
   real, dimension(:), pointer :: minimum
   integer :: ndims2, npoints2

   integer :: test, i, j

   do test = 1, ntests
      call random_number(r_ndims)
      call random_number(r_npoints)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))
      allocate(minimum(ndims))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange
      write(0,*) "Test number: ", test
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      tree => kd_construct(arry1(:,:))

      ! Test each dimension
      do test_dim = 1, ndims, 1
         call kd_find_min(tree, minimum, test_dim)

         if( .NOT. ( minimum(test_dim) == minval(arry1(test_dim,:)) )) then
            write(0,*) "We found the wrong minimum in dimension: ", test_dim
            write(0,*) "Min we found was: ", minimum
            write(0,*) "Actual minimum was: ", minval(arry1(test_dim,:))
            stop
         endif
      enddo

      call kd_free(tree)

      deallocate(arry1)
      deallocate(arry2)
      deallocate(minimum)
   enddo

   write(0,*) "All tests finished!"

end subroutine test_min

recursive subroutine print_tree(tree, left_or_right, level)

   implicit none
   type(kdnode), intent(inout), pointer :: tree
   integer, value :: left_or_right
   integer, value :: level

   if ( .NOT. associated(tree)) then
      write(0,*) " I AM AN UNASSOCIATED NODE"
      return
   endif

   if (associated(tree) .AND. associated(tree % point)) then
      if (left_or_right == -1) then
         write(0,*) "Level: ", level, "Left: ", tree % point(:) 
      elseif (left_or_right == 0) then
         write(0,*) "Level: ", level, "Root: ", tree % point(:)
      elseif (left_or_right == 1) then
         write(0,*) "Level: ", level, "Right: ", tree % point(:)
      endif
   endif

   if ( associated(tree % left)) call print_tree(tree % left, -1, level + 1)
   if ( associated(tree % right)) call print_tree(tree % right, 1, level + 1)

end subroutine print_tree

subroutine test_remove(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(ndims) :: min_point 
   real, dimension(3,10) :: arry3
   integer :: i, j
   real :: r_ndims, r_npoints
   real, dimension(ndims) :: minimum
   integer :: ndims2, npoints2

   do j = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", j
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange
      tree => kd_construct(arry1(:,:))

      do i = 1, size(arry1, dim=2), 1
         if(associated(tree)) then
            if(.NOT. kd_remove(tree, arry1(:,i))) then
               write(0,*) "We were not able to remove the point: ", arry1(:,i)
               stop
            endif
         else
            write(0,*) "The tree is no long associated"
         endif
      enddo

      if ( associated(tree) ) then
         write(0,*) "This tree was still associated :("
         write(0,*) tree % point
         stop
      endif

      deallocate(arry1)
      deallocate(arry2)
      call kd_free(tree)

   enddo

end subroutine test_remove

end program kd_tester

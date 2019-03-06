module kd_tree
    
   !***********************************************************************
   !
   !  module kd_tree
   !
   !>
   !
   !-----------------------------------------------------------------------
   implicit none

   private

   public :: kdnode

   ! Public Subroutines
   public :: kd_insert
   public :: kd_construct
   public :: kd_search
   public :: kd_remove
   public :: kd_free
   public :: kd_find_min

   type kdnode
      type (kdnode), pointer :: left => null()
      type (kdnode), pointer :: right => null()

      integer :: split_dim
      real, dimension(:), pointer :: point
   end type kdnode

   contains

   !***********************************************************************
   !
   !  routine kd_insert
   !
   !> This routine adds a point, `val(:)` to an existing KD-Tree or it creates 
   !> a new kdtree if kdtree is not associated.
   !>
   !> The `dim` variable is an optional dummy argument and should not be set 
   !> to zero, but preferably not used when starting with the root of a kdtree
   !> i.e.:
   !> 
   !> call kdtree(kdtree, val)
   !>
   !
   !-----------------------------------------------------------------------
   recursive subroutine kd_insert(kdtree, val, dim)

      implicit none

      ! Input Variables
      type(kdnode), intent(inout), pointer :: kdtree
      real, dimension(:), intent(in) :: val
      integer, optional, value :: dim

      integer :: d 

      if (.NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = modulo(d, size(val)) + 1

      if ( .NOT. associated(kdtree)) then
         allocate(kdtree)
         allocate(kdtree % point(size(val)))
         kdtree % left => null()
         kdtree % right => null()
         kdtree % split_dim = d
         kdtree % point(:) = val(:)
         return
      endif

      if (val(kdtree % split_dim ) > kdtree % point(kdtree % split_dim)) then
         call kd_insert(kdtree % right, val, kdtree % split_dim)
      else
         call kd_insert(kdtree % left, val, kdtree % split_dim)
      endif

   end subroutine kd_insert

   !***********************************************************************
   !
   !  recusrive routine kd_construct_internal
   !
   !> Recursive function for kd_construct. See kd_construct for
   !> more information.
   !>
   !
   !-----------------------------------------------------------------------
   recursive function kd_construct_internal(points, ndims, npoints, dim) result(tree)

      implicit none

      ! Input Varaibles
      real, dimension(:,:) :: points
      integer, intent(in) :: ndims
      integer, value :: npoints
      integer, value :: dim

      ! Return Value
      type (kdnode), pointer :: tree

      ! Local Variables
      integer :: median

      if (npoints < 1) then
         tree => null()
         return
      endif

      ! Sort the points at the split dimension
      dim = mod(dim, ndims) + 1
      call quickSort(points, dim, 1, npoints)

      median = (1 + npoints) / 2 

      allocate(tree) ! Allocate the node
      allocate(tree % point(ndims)) ! Allocate the point for that node
      tree % split_dim = dim
      tree % point = points(:,median)

      ! Build the right and left sub-trees but do not include the 
      ! node that was just allocated (i.e. points(:, median))
      tree % left => kd_construct_internal(points(:,1:median-1), ndims, median - 1, tree % split_dim)
      tree % right => kd_construct_internal(points(:,median+1:npoints), ndims, npoints - median, tree % split_dim)

   end function kd_construct_internal

   !***********************************************************************
   !
   !  routine kd_construct
   !
   !> This routine creates a balanced KD-Tree from a set of K-dimensional 
   !> points via quicksort and it returns a pointer to the root node of that
   !> tree. The points dummy argument, should be an array with the dimensions
   !> defined as: `points(k, n)` with k being the number of dimensions, and n
   !> being the number of points.    
   !>
   !> tree => kd_construct(points)
   !>
   !
   !-----------------------------------------------------------------------
   function kd_construct(points) result(tree)

      implicit none

      ! Input Varaibles
      real, dimension(:,:) :: points

      ! Return Value
      type (kdnode), pointer :: tree

      ! Local Varaibles
      integer :: ndims
      integer :: npoints

      ndims = size(points, dim=1)
      npoints = size(points, dim=2)
      
      if(npoints < 1) then
         ! No points were passed in, return null
         write(0,*) "ERROR: kd_tree - No points were passed in to construct!"
         tree => null()
         return
      endif

      tree => kd_construct_internal(points(:,:), ndims, npoints, 0)

   end function kd_construct


   !***********************************************************************
   !
   !  recursive routine kd_search_internal
   !
   !> Recursive subroutine for kd_search. See kd_search for more
   !> information.
   !
   !-----------------------------------------------------------------------
   recursive subroutine kd_search_internal(kdtree, query, res, distance)

      implicit none

      ! Input Variables
      type(kdnode), pointer, intent(in) :: kdtree
      real, dimension(:), intent(in) :: query
      real, dimension(:), intent(inout) :: res
      real, intent(inout) :: distance

      ! Local Values
      real :: current_distance

      current_distance = sum((kdtree % point(:) - query(:))**2)
      if (current_distance < distance) then
         distance = current_distance
         res = kdtree % point(:) 
      endif

      !
      ! To find the nearest point, we first attempt to find the point in the same manner
      ! as a single deminsion BST.
      !
      ! However, because we are looking for the nearest neighbor, then there might be
      ! a possibility that the nearest neighbor is on the otherside of the tree.
      !
      ! Thus, to determine if we need to search the opposite child we just searched, we
      ! will compare the distance of the current minimum distance, and the root node
      ! that we branched off of.
      !
      ! If the distance to the root node, is less then the current minimum distance,
      ! then the nearist neighbor might be in opposite child.
      !

      ! TODO: Double precision calculations

      if (query(kdtree % split_dim) > kdtree % point(kdtree % split_dim)) then
         if (associated(kdtree % right)) then ! Search right
            call kd_search_internal(kdtree % right, query, res, distance)
         endif
         if ((kdtree % point(kdtree % split_dim) - query(kdtree % split_dim))**2 <= distance .AND. associated(kdtree % left)) then 
            call kd_search_internal(kdtree % left, query, res, distance)
         endif
      else if (query(kdtree % split_dim) < kdtree % point(kdtree % split_dim)) then 
         if (associated(kdtree % left)) then ! Search left
            call kd_search_internal(kdtree % left, query, res, distance)
         endif
         if ((kdtree % point(kdtree % split_dim) - query(kdtree % split_dim))**2 <= distance .AND. associated(kdtree % right)) then
            call kd_search_internal(kdtree % right, query, res, distance)
         endif
      else ! Nearest point could be in either left or right subtree, so search both
         if(associated(kdtree % right)) call kd_search_internal(kdtree % right, query, res, distance)
         if(associated(kdtree % left)) call kd_search_internal(kdtree % left, query, res, distance)
      endif

   end subroutine kd_search_internal

   !***********************************************************************
   !
   !  routine kd_search
   !
   !> Find `point` within `kdtree` and return the nearest neighbor (or the point)
   !> within `result` return the distance in `min_d`.
   !>
   !
   !-----------------------------------------------------------------------
   subroutine kd_search(kdtree, query, res, distance)

      implicit none
      type(kdnode), pointer, intent(in) :: kdtree
      real, dimension(:), intent(in) :: query
      real, dimension(:), intent(inout) :: res
      real, intent(inout) :: distance

      if (size(kdtree % point) /= size(query)) then
         write(0,*) "ERROR: Searching a ", size(kdtree % point), "dimensional kdtree for a point that only"
         write(0,*) "ERROR: ", size(query), " dimensions. Please supply a point of equal"
         write(0,*) "ERROR: dimensions!"
         return
      endif

      call kd_search_internal(kdtree, query, res, distance)

   end subroutine kd_search

   !***********************************************************************
   !
   !  routine kd_find_min_internal
   !
   !> Find the minmum value that lies within the dimension `dim` and return the entire
   !> point within `point`. On first call, `minimum(:)`, should be set to huge(minimum(:)).
   !>
   !>
   !> NOTE: This is a internal function, and should not be called, use kd_find_min
   !> instead.
   !
   !-----------------------------------------------------------------------
   recursive subroutine kd_find_min_internal(kdtree, point, dim, minimum, depth)

      implicit none

      ! Input variables
      type (kdnode), pointer :: kdtree
      real, dimension(:), intent(inout) :: point
      integer, intent(in) :: dim
      real, dimension(:), intent(inout) :: minimum 
      integer, optional, value :: depth

      ! Local variables
      integer :: ndims
      integer :: d

      ndims = size(point)
      

      if ( .NOT. present(depth)) then
         d = 0
      else
         d = depth
      endif

      d = mod(d, ndims) + 1

      ! Base Case
      if ( .NOT. associated(kdtree)) then
         return
      endif


      if(kdtree % point(dim) < minimum(dim)) then
         minimum(:) = kdtree % point(:)
         point(:) = kdtree % point(:)
      endif

      ! If the current split dimension (d) is equal to the dimension we asked for (dim)
      ! then we know that the smallest point in this dimension is to the left
      !
      ! If the current split dimension (d) is not equal to the dimension we're searching
      ! on (dim), then the minimum can either be to the left, or the right.

      if (d == dim) then
         if (associated(kdtree%left)) then
            call kd_find_min_internal(kdtree%left, point, dim, minimum, d)
            return
         endif
      else
         call kd_find_min_internal(kdtree%left, point, dim, minimum, d)
         call kd_find_min_internal(kdtree%right, point, dim, minimum, d)
      endif

      ! Else we do not know here the smallest value lies, so we must recursivly search both
      ! subtrees.

   end subroutine kd_find_min_internal

   !***********************************************************************
   !
   !  routine kd_find_min
   !
   !> Find the minimum point within a KD-tree in the dimension `dim` and return
   !> that coordinate within `point`.
   !>
   !> If calling from a kd_tree routine, depth should be the split point of kdtree the
   !> tree minus 1. For Example:
   !>
   !> K = 3, Current Dimension = 1
   !>
   !> call kd_find(kdtree, point, dim, 0)
   !>
   !> or
   !>
   !> K = 3, Current Dimension = 1
   !> call kd_find(kdtree % right, point, dim, 1)
   !>
   !> Because kdtree % right is split upon along the 2nd dimension
   !>
   !
   !-----------------------------------------------------------------------
   subroutine kd_find_min(kdtree, point, dim, depth)

      implicit none

      ! Input Variables
      type (kdnode), pointer :: kdtree
      real, dimension(:), intent(inout) :: point
      integer, intent(in), value :: dim
      integer, optional, value :: depth

      ! Local variables
      real, dimension(size(point)) :: minimum
      integer :: ndims

      minimum = huge(minimum)
      ndims = size(kdtree % point, dim=1)

      if (ndims /= size(point)) then
         write(0,*) "ERROR: The kd Tree has ", ndims, " dimensions and the "
         write(0,*) "ERROR: variable to hold the return point (point) has "
         write(0,*) "ERROR: ", size(point), "dimensions."
         write(0,*) "ERROR: Please insure that they have the same dimensions"
         return
      endif

      if (.not. present(depth)) then
         call kd_find_min_internal(kdtree, point, dim, minimum)
      else
         call kd_find_min_internal(kdtree, point, dim, minimum, depth)
      endif

   end subroutine kd_find_min

   !***********************************************************************
   !
   !  routine kd_remove
   !
   !> Remove `point` from `kdtree` and return .TRUE. if the point was removed
   !> succesfully and .FALSE. if the point was unable to be removed.
   !>
   !> A point that is not a leaf node, will be removed and replaced with the 
   !> point that contains the lowest value of the deleted node's split dimension
   !> from the right subtree. If the removed point does not contain a right subtree
   !> then the minimum point will be replaced with the minimum found within the
   !> left subtree. The replacement node will then be recusrively deleted.
   !>
   !
   !-----------------------------------------------------------------------
   recursive function kd_remove(kdtree, point, dim) result(ierr)

      implicit none

      ! Input Variables
      type (kdnode), pointer :: kdtree
      real, dimension(:), intent(in) :: point
      integer, optional, value :: dim

      ! Return Value
      real, dimension(size(point)) :: min_point
      integer :: d
      integer :: ndims
      logical :: ierr

      ndims = size(point)

      if ( .NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = mod(d, ndims) + 1

      ierr = .FALSE.
         
      if ( .NOT. associated(kdtree)) then
         return
      endif

      if (ndims /= size(kdtree % point, dim=1)) then
         write(0,*) "ERROR: The KD Tree has ", ndims, " dimensions and the "
         write(0,*) "ERROR: requested point to remove has ", size(point), " dimensions"
         write(0,*) "ERROR: Please insure that they have the same dimensions"
         return
      endif

      ! If the point that was asked to remove is a leaf node, then delete it trivially
      if( .not. associated(kdtree % left) .AND. .not. associated(kdtree % right)) then
         if (all(kdtree % point(:) == point(:))) then
            deallocate(kdtree % point)
            deallocate(kdtree)
            ierr = .TRUE.
            return
         endif
      endif

      if(all(kdtree % point(:) == point(:))) then ! The current node equals the requested point

         !
         ! If the requested node is not a leaf node, then a replacement will
         ! need to be found for it.
         !
         ! Replace the deleted node with the smallest value of the current split
         ! dimension from the right subtree if possible and then recusrivly call
         ! remove on the replacment.
         !
         ! If there is no right subtree, then replace the deleted node with the
         ! point that contains the smallest value of the current split dimension
         ! and then swap the right and left subtrees. (i.e. the left becomes
         ! null). And similarly, call call remove on the replacement.
         !

         if (associated(kdtree % right)) then ! Find the minimum on the right side
            call kd_find_min(kdtree % right, min_point, d, d)
            kdtree % point(:) = min_point(:) 

            ! Call delete on the node that we found to be the minimum
            if ( .NOT. kd_remove(kdtree % right, min_point, d)) then
               ! TODO: Probably change this stop to something else
               stop
            endif
            ierr = .TRUE.
            return

         elseif ( associated(kdtree % left) ) then ! Find the minimum on the left side
            call kd_find_min(kdtree % left, min_point, d, d)
            kdtree % point(:) = min_point(:) 

            ! Call delete on the node that we found to be the minimum
            if ( .NOT. kd_remove(kdtree % left, min_point, d)) then
               ! TODO: Probably change this stop to something else
               stop
            endif

            ! If we replaced with the left subtree, then 
            kdtree % right => kdtree % left
            kdtree % left => null()
            ierr = .TRUE.
            return
         endif
      endif

      ! If kdtree % point(:) is not the current node, then we still need to search
      ! for it. So make a deicion based on the current split dimension

      ! Search Left
      if (point(d) < kdtree % point(d) .AND. associated(kdtree % left)) then
         if (kd_remove(kdtree % left, point, d)) then
            ierr = .TRUE.
            return
         endif
      endif
      ! Search Right
      if (point(d) > kdtree % point(d) .AND. associated(kdtree % right)) then
         if (kd_remove(kdtree % right, point, d)) then
            ierr = .TRUE.
            return
         endif
      endif

      ! If we fall off the end, then the point request is mot likely not within the
      ! tree

   end function kd_remove

   !***********************************************************************
   !
   !  routine kd_free
   !
   !> Recursivly deallocate all nodes within `kdtree` including `kdtree` itself.
   !>
   !
   !-----------------------------------------------------------------------
   recursive subroutine kd_free(kdtree)

      implicit none
      type(kdnode), pointer :: kdtree

      if (.not. associated(kdtree)) then
         return
      endif

      if (associated(kdtree % left)) then
         call kd_free(kdtree % left)
      endif

      if (associated(kdtree % right)) then
         call kd_free(kdtree % right)
      endif

      deallocate(kdtree % point)
      deallocate(kdtree)

   end subroutine kd_free


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sorts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !***********************************************************************
   !
   !  routine kd_quicksort
   !
   !> Sort points starting from arrayStart, to arrayEnd along the given dimension 
   !> `dim`. If two points are swapped, the entire K-Coordinate point are swapped.
   !>
   !> TODO: Change the name of this function to kd_quicksort
   !
   !-----------------------------------------------------------------------
   !recursive subroutine kd_quicksort(array, dim, arrayStart, arrayEnd)
   recursive subroutine quickSort(array, dim, arrayStart, arrayEnd)

      implicit none

      ! Input Variables
      real, dimension(:,:) :: array
      integer, intent(in), value :: dim
      integer, intent(in), value :: arrayStart, arrayEnd

      ! Local Variables
      integer :: ndims, npoints
      real, dimension(size(array, dim=1)) :: temp
      real, dimension(size(array, dim=1)) :: pivot_value

      integer :: l, r, pivot, s

      ndims = size(array, dim=1)
      npoints = arrayEnd

      if ((arrayEnd - arrayStart) < 1) then
         return
      endif

      ! Create the left, right, and start pointers
      l = arrayStart
      r = arrayEnd - 1
      s = l

      pivot = (l+r)/2
      pivot_value = array(:, pivot)

      ! Move the pivot to the far right
      temp(:) = array(:,pivot)
      array(:,pivot) = array(:,arrayEnd)
      array(:,arrayEnd) = temp(:)

      do while ( .TRUE. )
         ! Advance the left pointer until it is a value less then our pivot_value(dim)
         do while ( .TRUE. )
            if (array(dim, l) < pivot_value(dim)) then
               l = l + 1
            else
               exit
            endif
         enddo

         ! Advance the right pointer until it is a value more then our pivot_value(dim)
         do while ( .TRUE. )
            if ( r <= 0 ) then
               exit
            endif

            if(array(dim, r) > pivot_value(dim)) then
               r = r - 1
            else
               exit
            endif
         enddo

         if ( l >= r ) then 
            exit
         else ! Swap elements about the pivot
            temp = array(:,l)
            array(:,l) = array(:,r)
            array(:,r) = temp
         endif
      enddo

      ! Move the pivot to l ended up
      temp(:) = array(:,l)
      array(:,l) = array(:,arrayEnd)
      array(:,arrayEnd) = temp(:)

      !Quick Sort on the lower partition
      call quickSort(array(:,:), dim, s, l-1)

      !Quick sort on the upper partition
      call quickSort(array(:,:), dim, l+1, arrayEnd)

   end subroutine quicksort

end module kd_tree

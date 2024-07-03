module vector_operations
   implicit none
   private
   public :: dotProduct, crossProduct, normTwo, dirCos

contains

   ! Calculates the dot product of two vectors
   function dotProduct(vec1, vec2) result(dot)
      real, intent(in), dimension(:) :: vec1, vec2
      real :: dot
      dot = sum(vec1 * vec2)
   end function dotProduct

   ! Computes the cross product of two three-dimensional vectors
   function crossProduct(vec1, vec2) result(cross)
      real, intent(in), dimension(3) :: vec1, vec2
      real, dimension(3) :: cross
      cross(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
      cross(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
      cross(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
   end function crossProduct

   ! Calculates the norm (2-norm) of a vector
   function normTwo(vec) result(norm)
      real, intent(in), dimension(:) :: vec
      real :: norm
      norm = sqrt(dotProduct(vec, vec))
   end function normTwo

   ! Computes the direction cosines of a vector
   function dirCos(vec) result(cosines)
      real, intent(in), dimension(3) :: vec
      real, dimension(3) :: cosines
      real :: vecNorm
      vecNorm = normTwo(vec)
      if (vecNorm > 0.0) then
         cosines = vec / vecNorm
      else
         cosines = 0.0
      end if
   end function dirCos

end module vector_operations

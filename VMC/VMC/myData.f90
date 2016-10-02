module myData
    type line_info
        real :: time
        real :: temp
    end type

    type datalink
        type(line_info) :: time_temp
        type(datalink), pointer :: next
    end type

    type(line_info), allocatable :: FileInfoArr(:)
end module
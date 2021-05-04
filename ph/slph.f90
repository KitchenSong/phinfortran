program slph

use util
use config
use view_struct
use dyn

implicit none

call load_configure()
call gen_struct() 
call gen_dyn() 

end program

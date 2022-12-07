#!/bin/bash
touch par_for_output_$1.log
../CoalBathes par CoalData_$1.root >> par_for_output_$1.log 
echo "Done"

#!/bin/bash
touch par_for_output_$1.log
../CoalBathes par_$1 CoalData_$1.root >> par_for_output_$1.log 
echo "Done"

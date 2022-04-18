echo '<Info>'
echo 'Co, Cf = Initial and Final Cost'
echo 'Wo, Wf = Initial and Final Worst Case'
echo 'dC, dW = Cf-Co and Wf-Wo'
echo 'nS = Number of Swaps'
echo 'aC = All Costs'
echo ''

for i in {1..1000}; do

	python3 src/main.py src/json2graph.py src/mappingGRN.py
	echo ''

done

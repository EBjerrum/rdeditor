from mendeleev import element


elements={}
symboltoint = {}


for i in range(1,119):
	e = element(i)
	elements[i] = {"Name":e.name, "Symbol":e.symbol,"Group": e.group_id,"Period":e.period}
	symboltoint[e.symbol] = i




#http://www.sfu.ca/%7Evdabbagh/gap/doc/chap0.html
LoadPackage("repsn");;



#Matrix forms of the gens of BO in our chosen rep
IdMat := [
           [1,0],
           [0,1]
          ];;


negIdMat := [
           [-1,0],
           [0,-1]
             ];;


iMat := [
          [0,E(4)],
          [E(4),0]
         ];;


jMat := [
          [0,1],
          [-1,0]
         ];;


kMat := [
          [E(4),0],
          [0,-E(4)]
         ];;


uMat := -1/2 * (IdMat + iMat + jMat + kMat);;


tMat := 1/Sqrt(2) * (IdMat + iMat);;



#Our 'canonical' matrix form of BO, which happens to be the fundemental rep.
#We only need this to test that we got the right elements from GAP's basis for BO.
canBO := Group([negIdMat, jMat, kMat, uMat, tMat], IdMat);;

#GAP's internal form of BO
BO := SmallGroup(48,28);;


#From repsn package, Irr sends a group to its character table in the form of a list.
BOCharTable := Irr(BO);;

#We split out all the representations BO -> V
irr1 := IrreducibleAffordingRepresentation(BOCharTable[1]);;
irr2 := IrreducibleAffordingRepresentation(BOCharTable[2]);;
irr3 := IrreducibleAffordingRepresentation(BOCharTable[3]);;
irr4 := IrreducibleAffordingRepresentation(BOCharTable[4]);;
irr5 := IrreducibleAffordingRepresentation(BOCharTable[5]);;
irr6 := IrreducibleAffordingRepresentation(BOCharTable[6]);;
irr7 := IrreducibleAffordingRepresentation(BOCharTable[7]);;
irr8 := IrreducibleAffordingRepresentation(BOCharTable[8]);;

#Here we single out the fundemental rep, and pick out the elements we used to parameterize BO.
BOFunRep := List(BO, x -> rec(g := x, gmat := x^irr4));;
id := Filtered(BOFunRep, x->x.gmat = IdMat)[1].g;;
nid := Filtered(BOFunRep, x->x.gmat = negIdMat)[1].g;;
i := Filtered(BOFunRep, x->x.gmat = iMat)[1].g;;
j := Filtered(BOFunRep, x->x.gmat = jMat)[1].g;;
k := Filtered(BOFunRep, x->x.gmat = kMat)[1].g;;
u := Filtered(BOFunRep, x->x.gmat = uMat)[1].g;;
t := Filtered(BOFunRep, x->x.gmat = tMat)[1].g;;



# Generate all of the bitstrings corresponding to valid elements.
bitStrings := Cartesian(ListWithIdenticalEntries(6,[0,1]));;
bitStrings := Filtered(bitStrings, xs -> not(xs[4] = 1 and xs[5] = 1));


groupTable := [];;

for bitString in bitStrings do
    b1 := bitString[1];;
    b2 := bitString[2];;
    b3 := bitString[3];;
    b4 := bitString[4];;
    b5 := bitString[5];;
    b6 := bitString[6];;

    #Used for testing against irr4
    #canEl := negIdMat^b1 * jMat^b2 * kMat^b3 * uMat^(2*b4 + b5) * tMat^b6;;

    el := nid^b1 * j^b2 * k^b3 * u^(2*b4 + b5) * t^b6;;
    el1 := el^irr1;;
    el2 := el^irr2;;
    el3 := el^irr3;;
    el4 := el^irr4;;
    el5 := el^irr5;;
    el6 := el^irr6;;
    el7 := el^irr7;;
    el8 := el^irr8;;

    conjClass := ConjugacyClass(BO,el);;
    classRep := Representative(conjClass);;
    tr1 := Trace(classRep^irr1);;
    tr2 := Trace(classRep^irr2);;
    tr3 := Trace(classRep^irr3);;
    tr4 := Trace(classRep^irr4);;
    tr5 := Trace(classRep^irr5);;
    tr6 := Trace(classRep^irr6);;
    tr7 := Trace(classRep^irr7);;
    tr8 := Trace(classRep^irr8);;

    Add(groupTable, rec(
                         address := bitString,
                         intAddress := b6 + 2*b5 + 4*b4 + 8*b3 + 16*b2 + 32*b1,
                         #canEl := canEl,
                         el := el,
                         1 := el1,
                         2 := el2,
                         3 := el3,
                         4 := el4,
                         5 := el5,
                         6 := el6,
                         7 := el7,
                         8 := el8,
                         class := conjClass,
                         tr1 := tr1,
                         tr2 := tr2,
                         tr3 := tr3,
                         tr4 := tr4,
                         tr5 := tr5,
                         tr6 := tr6,
                         tr7 := tr7,
                         tr8 := tr8,

                        ));;
od;;

#Agrees with our definition of the generators?
#To use this to verify, uncomment the lines canEl := in the bitString loop above
#test := ForAll(groupTable, x -> x.4 = x.canEl);;


# From our table we can find which elements have are in each class in our character table
classes := ConjugacyClasses(BO);;
class1 := Filtered(groupTable, x -> x.class = classes[1]);;
class2 := Filtered(groupTable, x -> x.class = classes[2]);;
class3 := Filtered(groupTable, x -> x.class = classes[3]);;
class4 := Filtered(groupTable, x -> x.class = classes[4]);;
class5 := Filtered(groupTable, x -> x.class = classes[5]);;
class6 := Filtered(groupTable, x -> x.class = classes[6]);;
class7 := Filtered(groupTable, x -> x.class = classes[7]);;
class8 := Filtered(groupTable, x -> x.class = classes[8]);;

#We can extract specifically list of addresses of elements in each class
elInClass1 := List(class1, x -> x.intAddress);;
elInClass2 := List(class2, x -> x.intAddress);;
elInClass3 := List(class3, x -> x.intAddress);;
elInClass4 := List(class4, x -> x.intAddress);;
elInClass5 := List(class5, x -> x.intAddress);;
elInClass6 := List(class6, x -> x.intAddress);;
elInClass7 := List(class7, x -> x.intAddress);;
elInClass8 := List(class8, x -> x.intAddress);;
Sort(elInClass1);;
Sort(elInClass2);;
Sort(elInClass3);;
Sort(elInClass4);;
Sort(elInClass5);;
Sort(elInClass6);;
Sort(elInClass7);;
Sort(elInClass8);;

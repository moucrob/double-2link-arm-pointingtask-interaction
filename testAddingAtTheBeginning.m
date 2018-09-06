function vecp1 = testAddingAtTheBeginning(elementToAdd,vec)

tmp = fliplr(vec);
tmp(end+1) = elementToAdd;
vecp1 = fliplr(vec);

end
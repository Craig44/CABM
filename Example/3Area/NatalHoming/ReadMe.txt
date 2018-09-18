Some issues I have come across with implementing Natal Homeing, CASAL does not allow an area
to be a source cell for a single stock in a migration process. That is you cannot simultaneously define
the movement from EN -> HG aswell as EN -> BP so what they do in CASAL is define a two step movement
EN -> HG and then HG -> BP. For cerain combinations this actually seems sensible (in a biological sense) 
but for others such as the source cell = HG it doesn't for example you would have to describe the movement
from HG -> EN and then EN -> BP which doesn't make sense at all. Moral of the story is I don't like this, and
will only do it to test code otherwise we will not repeat this hack.
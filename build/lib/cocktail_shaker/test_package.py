# Implemented Strictly for Testing
# --------------------------------
if __name__ == '__main__':

    from cocktail_shaker.file_handler import FileWriter
    from cocktail_shaker.functional_group_enumerator import Cocktail

    cocktail = Cocktail(['c1cc(CCCO)ccc1'])
    compounds_result = cocktail.shake(functional_groups=['Azides'])

    FileWriter('testst', compounds_result, 'sdf', fragementation=)

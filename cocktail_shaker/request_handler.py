#!/usr/bin/env python
#
# Request Handling for MolPort and Cactus
#
# ----------------------------------------------------------

class RequestError(Exception):

    __version_error_parser__ = "1.1.0"
    __allow_update__ = False

    """

    Raise a Request Handler Error if the response code is not 200. 

    """
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors

class CactusRequestHandler(object):


    __version__ = "1.1.0"
    __allow_update__ = False

    """
    
    Handler API requests to external libraries
    
    
    """

    def __init__(self, url):

        self.url = url

    def get(self):
        """

        Issues Requests specifically for the Cactus website.

        """

        # imports
        # -------
        from urllib.request import urlopen

        try:
            request = urlopen(self.url)
        except RequestError as e:
            print ("StackTrace: %s" % e)
            print ("Error handling request, please try again or contact Suliman Sharif")

        if request.getcode() != 200:
            raise RequestError(message="Error handling request, please try again or contact Suliman Sharif",
                               errors='Request Error')

        return request

class Resolver(object):

    __version__ = "1.1.0"
    __allow_update__ = False

    """
    
    Handle responses coming back from a varied type of requests dependent. 
    
    """

    def __init__(self, request):

        self.request = request

    def cactus_mol2_resolver(self):

        from lxml import etree
        response = etree.parse(self.request).getroot()
        for child in response.iter('item'):
            response = child.text

        return response





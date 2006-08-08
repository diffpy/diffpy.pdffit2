/***********************************************************************
* Short Title: unit tests driver
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace std;
using namespace CppUnit;

int main()
{
    // Get the top level suite from the registry
    Test *suite = TestFactoryRegistry::getRegistry().makeTest();

    // Adds the test to the list of test to run
    TextUi::TestRunner runner;
    runner.addTest( suite );

    // Change the default outputter to a compiler error format outputter
    runner.setOutputter( new CompilerOutputter(&runner.result(), cerr) );
    // Run the tests.
    bool wasSucessful = runner.run();

    // Return error code 1 if the one of test failed.
    return wasSucessful ? EXIT_SUCCESS : EXIT_FAILURE;
}

/***********************************************************************
* Too see what people have been up to just run:
*   cvs log
***********************************************************************/

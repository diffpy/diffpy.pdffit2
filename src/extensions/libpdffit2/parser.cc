/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2006 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Jacques Bloch
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* Formula parser implemented as methods of Fit class
*
* Comments:
*
***********************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

#include "MathUtils.h"
#include "pdffit.h"

static vector<double> vstack;
static vector<vector<double> > dstack;   // stack with derivative vectors
static bool deriv;

string stackvar(int i)
{
    ostringstream expr;
    expr << '#' << i;
    return(expr.str());
}

/*****************************************************************************************
  Searches for all parameters in the expression and replace them by stack pointers
 *****************************************************************************************/
string Fit::substitute_pars(string &expr)
{
    int ipos=-1, bopen, bclose;
    unsigned int i;
    string idexpr;

    // search the expression for the next @-sign
    while((bopen=expr.find('@',ipos+1)) != (int) string::npos)
    {
	int id;

	istringstream inexpr(expr);
	inexpr.seekg(bopen+1);

	inexpr >> id;

	if (!inexpr)
	    throw parseError(expr);

	bclose = inexpr.tellg();

	// expression must be a valid parameter id
	int ipar = parfind(id);

	if (ipar == -1)
	{
	    ostringstream msg;
	    msg << "parameter " << id << " undefined";
	    throw constraintError(msg.str());
	}

	// check if same parameter already occurred in this formula
	for (i=0; i<used.size(); i++)
	{
	    if(used[i] == ipar)
		break;
	}
	if (i==used.size())
	{
	    used.push_back(ipar);
	    vstack.push_back(p[ipar]);
	}

	expr.replace(bopen,bclose-bopen,stackvar(i));

    }

    // check if any derivatives are required
    deriv = false;
    for (i=0; i<used.size(); i++)
	deriv = deriv || ip[used[i]];

    // make a derivative array with a 1 corresponding to this parameter
    // only make derivatives if at least 1 parameter will be refined
    if (deriv)
    {
	dstack.clear();

	vector<double> dnumdp(used.size());
	for (i=0; i<used.size(); i++)
	{
	    if (ip[used[i]])
	    {
		dnumdp[i] = 1;
		dstack.push_back(dnumdp);
		dnumdp[i] = 0;
	    }
	    else
		// if parameter not refined: push empty derivative vector on stack
		dstack.push_back(vector<double>());
	}
    }

    return expr;  // okay
}

// get next number from expression
// only fill derivatives if that number has derivatives
double Fit::getnum(istringstream &inexpr, vector<double> &dnumdp)
{
    double num;
    char c;
    int id;
    map<string,Builtin>::iterator iter;

    // Various possibilities when reading the next value inside the
    // unbracketted expression

    if (deriv) dnumdp.clear();

    inexpr >> c;

    if (inexpr.eof())
	throw parseError("Error while reading value");

    // if first character is a minus sign, call getnum again with expression stripped
    // from its minus sign
    if ((c=='-') || (c=='+'))
    {
	num = getnum(inexpr, dnumdp);
	if (deriv && !dnumdp.empty() && (c=='-'))
	{
	    for (unsigned int i=0; i<used.size(); i++)
		dnumdp[i] = -dnumdp[i];
	}
	if (c=='-')
	    return -num;
	else
	    return num;
    }
    // if first character is a letter -> check for built-in function
    // read-up to next #
    else if ((c>='a') && (c<='z'))
    {
	int start=inexpr.tellg();
	start--;
        string::size_type end;
        end = inexpr.str().find('#', start);
	if (end == string::npos)
        {
	    throw parseError("Error while reading builtin function arguments");
        }

	string sbuiltin = inexpr.str().substr(start,end-start);
	if ((iter=builtin.find(sbuiltin)) == builtin.end())
        {
	    throw parseError("Unknown builtin function");
        }

	// set read pointer behind #-sign
	inexpr.seekg(end+1);

	// read argument number
	inexpr >> id;
	if (!inexpr)
	    throw parseError(inexpr.str());
	num = iter->second.func(vstack[id]);

	// fill the partial derivatives if function argument had derivatives
	if (deriv && !dstack[id].empty())
	{
	    double fder=iter->second.deriv(vstack[id]);
	    dnumdp = vector<double>(used.size());
	    for(unsigned int i=0; i<used.size(); i++)
		dnumdp[i] = fder*dstack[id][i];
	}
    }
    else if (c == '#')
    {
	inexpr >> id;
	if (!inexpr)
	    throw parseError(inexpr.str());
	num = vstack[id];
	if (deriv && !dstack[id].empty()) dnumdp = dstack[id];
    }
    else
    {
	inexpr.unget();
	inexpr >> num;

	if (!inexpr)
	    throw parseError("Error while reading value");
    }
    return num;
}

/**********************************************
  computes a basic expression without brackets
 ***********************************************/
double Fit::compute(string &expr, vector<double> &dnumdp)
{
    ostringstream ostreamexpr;
    string opstring="*/+-";
    double num1, num2, num;
    vector<double> dnum1dp, dnum2dp;
    char op;  // operator

    // transform string to string stream to read numbers more easily
    istringstream inexpr(expr);

    // Two parsing passes are performed to ensure correct precedence of
    // operators: first for * and /, then for + and -
    int pass = 1;

    while(1)
    {
	int pos, end;

	if (deriv) dnumdp.clear();

	pos = inexpr.tellg();

	num1 = getnum(inexpr, dnum1dp);

	inexpr >> op;

	// check if this ends the expression
	// if pass==2: expression completely parsed
	// if pass==1: go to pass 2
	if (inexpr.eof())
	{
	    if (pass==2)
	    {
		num = num1;
		if (deriv && !dnum1dp.empty())
		    dnumdp = dnum1dp;
		break;
	    }
	    else
	    {
		pass = 2;
		inexpr.clear();
		inexpr.str(expr);
		continue;
	    }
	}

	if (!inexpr)
	    throw parseError(expr);

	if (opstring.find(op) == string::npos)
	    throw parseError(expr);

	// if a + or - operation is encountered in pass 1 -> look for next operation
	if ( (pass==1) && ((op == '+') || op == '-') )
	{
	    continue;
	}

	num2 = getnum(inexpr,dnum2dp);

	if (!dnum1dp.empty() || !dnum2dp.empty())
	    dnumdp = vector<double>(used.size());

	switch (op)
	{
	    case '+':
		num = num1+num2;

		if (deriv)
		{
		    if (!dnum1dp.empty()) dnumdp = dnum1dp;
		    if (!dnum2dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] += dnum2dp[i];
		}
		break;

	    case '-':
		num = num1-num2;

		if (deriv)
		{
		    if (!dnum1dp.empty()) dnumdp = dnum1dp;
		    if (!dnum2dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] -= dnum2dp[i];
		}
		break;

	    case '*':
		num = num1*num2;

		if (deriv)
		{
		    if (!dnum1dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] += dnum1dp[i]*num2;
		    if (!dnum2dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] += num1*dnum2dp[i];
		}
		break;

	    case '/':
		num = num1/num2;

		if (deriv)
		{
		    if (!dnum1dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] += dnum1dp[i]/num2;
		    if (!dnum2dp.empty())
			for (unsigned int i=0; i< used.size(); i++)
			    dnumdp[i] -= num1*dnum2dp[i]/sqr(num2);
		}
		break;
	}

	end = inexpr.tellg();

	// replace computed substring by stack pointer, and push value on the stack
	expr.replace(pos,end-pos,stackvar(vstack.size()));
	vstack.push_back(num);

	// store the derivatives as well. Note that the derivative vector may be empty
	//  if no derivatives was non-zero.
	if (deriv)
	    dstack.push_back(dnumdp);

	inexpr.clear();
	inexpr.str(expr);
	inexpr.seekg(pos);
    }

    // the expression has been reduced to just one number
    return num;
}

/***************************************************************************************
  Searches the formula for bracketted expressions replace them by their numerical values
 ****************************************************************************************/

double Fit::parse(string line, vector<double> &dnumdp)
{
    int bopen, bclose;
    string expression;
    double num;

    if (line.find('#') != string::npos)
	throw parseError("Illegal character in formula");

    // the first elements on the stack will be the used parameter values
    vstack.clear();
    used.clear();

    // Search for parameters, check their existence and put on the stack
    line = substitute_pars(line);

    // look for first closing bracket
    while( (bclose=line.find(')')) != (int) string::npos)
    {
	// look for closest opening bracket in front of closing bracket
	if ( (bopen=line.rfind('(',bclose-1)) == (int) string::npos)
	    throw parseError("Unmatched brackets in formula");

	// isolate expression within the brackets
	expression = line.substr(bopen+1,bclose-bopen-1);

	// compute isolated expression
	num = compute(expression, dnumdp);

	// put the bracket-value on the stack and replace substring by stack pointer,
	// and push new value on the stack
	line.replace(bopen,bclose-bopen+1,stackvar(vstack.size()));
	vstack.push_back(num);

	if (deriv)
	    dstack.push_back(dnumdp);

    }

    // After all brackets have been worked out, the final expression is evaluated
    num = compute(line, dnumdp);

    return num;  // okay
}

// End of file

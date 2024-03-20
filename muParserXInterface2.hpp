#ifndef HH_MUPARSERXINTERFACE_HH
#define HH_MUPARSERXINTERFACE_HH

#include "mpParser.h"
#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace MuParserInterface
{
    /*! An interface to MuParserX to define a function

    It defines a functor representing a function \f$ \mathbb{R}^N \f$ to \f$ \mathbb{R}^M\f$.
    The input variables are defined as x[0], x[1], etc., and one can use the MuParserX
    syntax to create the expression.

    I assume that at compile time we know the size of the argument of the function
    and I keep the return value as template parameter. By default, both input and
    output are std::vector<double>.
    The input variables are indicated by x[]. An example of a valid expression:
    {sin(x[0]), x[1]*x[2]}

    @tparam N The number of arguments of the function we want to represent
    @tparam M The number of elements in the output vector
    @tparam ArgumentType any type that supports the addressing ([]) operator
    */
    template <int N, int M, class ArgumentType = std::array<double, N>>
    class muParserXInterface
    {
    public:
        //! Default constructor
        //!
        //! mup::pckALL_NON_COMPLEX|mup::pckMATRIX means that I do not want the module
        //! for complex numbers but I want to treat arrays and matrices in MuParserX
        //! expressions
        muParserXInterface()
            : My_e(),
              M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX), M_values(N, std::vector<double>(M, 0.0))
        {
            M_parser.DefineVar("x", mup::Variable(&M_values[0][0]));
        }

        //! Constructor that takes a string containing MuParserX expression
        muParserXInterface(const std::string expression) : muParserXInterface()
        {
            My_e = expression;
            M_parser.SetExpr(My_e.c_str());
        }

        /*!
        * The copy constructor
        *
        * MuParserX has a particular design, which obliges to define a special copy
        * constructor. The reason is that a MuParser engine stores the address of the
        * variables. So a normal copy would do a shallow copy, which is NOT what you
        * want. Moreover, because of a poor design, you may lose the expression.
        * That's why I keep a copy in the class as a string and redefine it in the
        * MuParser engine.
        *
        * @param mpi the muParserXInterface to be copied
        */
        muParserXInterface(muParserXInterface const &mpi)
            : My_e(mpi.My_e),
              M_parser(mup::pckALL_NON_COMPLEX | mup::pckMATRIX), M_values(N, std::vector<double>(M, 0.0))
        {
            M_parser.DefineVar("x", mup::Variable(&M_values[0][0]));
            M_parser.SetExpr(My_e.c_str());
        }

        /*!
        * The copy assignment operator
        *
        * MuParserX has a particular design, which obliges to define a special copy
        * assignment operator.
        * @param mpi the muParserXInterface to be copied
        * The copy constructor
        */
        muParserXInterface
        operator=(muParserXInterface const &mpi)
        {
            if (this != &mpi)
            {
                this->My_e = mpi.My_e;
                this->M_parser.ClearVar(); // clear the variables!
                this->M_values = mpi.M_values;
                M_parser.DefineVar("x", mup::Variable(&M_values[0][0]));
                M_parser.SetExpr(My_e.c_str());
            }
            return *this;
        }

        //! Sets the MuParserX expression.
        /*!
        * Beware, the input variables are indicated by x[].
        * Example of a valid expression: {sin(x[0]), x[1]*x[2]}
        * @par e The expression
        */
        void set_expression(const std::string &e)
        {
            My_e = e;
            M_parser.SetExpr(e.c_str());
        }

        auto operator()(ArgumentType const &x) const
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < M; ++j)
                {
                    M_values[i][j] = x[i][j];
                }
            }

            std::vector<std::vector<double>> ans(N, std::vector<double>(M));
            try
            {
                for (int i = 0; i < N; ++i)
                {
                    M_parser.SetExpr(My_e.c_str());
                    ans[i][0] = M_parser.Eval().GetFloat();
                }
            }
            catch (mup::ParserError &error)
            {
                std::cerr << "MuParserX error with code:" << error.GetCode() << std::endl;
                std::cerr << "While processing expression: " << error.GetExpr() << std::endl;
                std::cerr << "Error Message: " << error.GetMsg() << std::endl;
                throw error;
            }
            return ans;
        }

    private:
        // A copy of the MuParserX expression, used for the copy operations
        std::string My_e;
        // The MuParserX engine
        mup::ParserX M_parser;
        // The MuParserX values used to set the variables in the engine
        mutable std::vector<std::vector<double>> M_values;
    };
} // namespace MuParserInterface

#endif

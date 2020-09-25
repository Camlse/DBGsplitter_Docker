/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef DBGSPLITTER_DIFFERENTIAL_ANALYSIS_H
#define DBGSPLITTER_DIFFERENTIAL_ANALYSIS_H

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
namespace differential_analysis {
    class differential_analysis : public Tool {
    private:
        void populateParser();
        int nbReplicates;

    public:
        // Constructor
        differential_analysis();

        // Actual job done by the tool is here
        void execute();

        //overriding this in order to exit the tool when finding a problem with the arguments
        IProperties *run(int argc, char *argv[]) {
          IProperties *toReturn = Tool::run(argc, argv);
          if (!toReturn)
            std::exit(1);
          return toReturn;
        }
    };
}



#endif //DBGSPLITTER_DIFFERENTIAL_ANALYSIS_H

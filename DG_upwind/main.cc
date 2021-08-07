#include "source/DG_upwind.cc"


int main(int argc, char** argv)
{
  try
    {
	  std::string par_name = "";
	  if (argc > 1){
	    par_name = argv[1];
	  }
//	  deallog.depth_console(2); //solver infos
      AdvectionProblem<2> dgmethod;
	  if (par_name!=""){
		 dgmethod.initialize(par_name);
	  }


      dgmethod.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;


}

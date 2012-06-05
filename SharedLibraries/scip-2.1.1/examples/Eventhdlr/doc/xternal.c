/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Event Handler Example
 * @author   Stefan Heinz
 *
 * This example illustrates the use of an event handler within <a href="http://scip.zib.de">SCIP</a>. It extends the
 * default plugins of <a href="http://scip.zib.de">SCIP</a> with two additional plugin. That are, an event handler which
 * reacts on new best solutions and an event handler that acts on processed nodes. You should also read the section <a
 * href="http://scip.zib.de/doc/html/EVENT.html">How to add event handler</a> in the <a
 * href="http://scip.zib.de/doc/html/index.html">SCIP doxygen</a> documentation which explains the event handling in
 * general.
 *
 * The event handler event_bestsol.c and event_boundwriting.c shows how the create a customized event handler in <a
 * href="http://scip.zib.de">SCIP</a>. In case of an event handler there are two importand questions to answer. First
 * when to install the event handler and second when to remove the event handler.
 *
 * <b>Note:</b> You can replace the event type in this example with any none variable event. See type_event.h for a
 * complete list.
 *
 * @section INSTALL Installing the event handler
 *
 * In the following we describe the implementation of the specific callback methods on the event_bestsol.c event
 * handler.
 *
 * In case of this example, where we want to install an event handler which reacts on new best solution, we install the
 * event handler in the callback SCIP_DECL_EVENTINIT. At that point the problem was tranformed. Note that there are
 * heuristics which are called before even the presolving starts. For example heur_trivial.c. This means the callback
 * SCIP_DECL_EVENTINTSOL is to late to install the solution event since at that point we might already missed some
 * solutions.
 *
 * \code 
 * static
 * SCIP_DECL_EVENTINIT(eventInitBestsol)
 * {
 *    assert(scip != NULL);
 *    assert(eventhdlr != NULL);
 *    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
 *
 *    SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );
 *
 *    return SCIP_OKAY;
 * }
 * \endcode 
 *
 * The method SCIPcatchEvent() notifies <a href="http://scip.zib.de">SCIP</a> that we want to react on the event type best solution found.
 *
 * @section REMOVE Remove the event handler
 *
 * With respect to dropping the event handling, we perform that action in SCIP_DECL_EVENTEXIT which is the counter part
 * of SCIP_DECL_EVENTINIT. The callback SCIP_DECL_EVENTEXITSOL which is the counter part to SCIP_DECL_EVENTINTSOL does
 * not work. Since this method is called before the branch-and-bound process is freed. This, however, is not only done
 * after the problem is solved. It also happens when a restart is performed. Dropping the event handler within this
 * method would lead to the fact that we miss all best solution detected after the first restart. Therefore, the right
 * place to drop that event handler is the callback SCIP_DECL_EVENTEXIT. Below you find the source code.
 *
 * \code 
 * static
 * SCIP_DECL_EVENTEXIT(eventExitBestsol)
 * { 
 *    assert(scip != NULL);
 *    assert(eventhdlr != NULL);
 *    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
 *   
 *    SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );
 *    return SCIP_OKAY;
 * } 
 * \endcode
 *
 * The method SCIPdropEvent() tells SCIP that we want to drop the event type SCIP_EVENTTYPE_BESTSOLFOUND of belonging
 * to the your event handler.
 *
 * @section REACT React on events
 *
 * In the callback SCIP_DECL_EVENTEXEC, which is the execution method of the event handler, you see how to access the
 * current best solution. 
 *
 * \code
 * static
 * SCIP_DECL_EVENTEXEC(eventExecBestsol)
 * {  
 *   SCIP_SOL* bestsol;
 *   SCIP_Real solvalue;
 *
 *   assert(eventhdlr != NULL);
 *   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
 *   assert(event != NULL);
 *   assert(scip != NULL);
 *   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);
 *
 *   SCIPdebugMessage("exec method of event handler for best solution found\n");
 *  
 *   bestsol = SCIPgetBestSol(scip);
 *   assert(bestsol != NULL);
 *   solvalue = SCIPgetSolOrigObj(scip, bestsol);
 *  
 *   SCIPinfoMessage(scip, NULL, "found new best solution with solution value <%g>\n", solvalue);
 *  
 *   return SCIP_OKAY;
 * }
 * \endcode
 *
 *
 * @section ADDING Including the event handler plugin
 * 
 * <a href="http://scip.zib.de">SCIP</a> is plug and play based. This means, all plugins which should be used have be
 * included into the <a href="http://scip.zib.de">SCIP</a> environment. In the case of the event handler, we are duing
 * this after the <a href="http://scip.zib.de">SCIP</a> environment was created (see cmain.c).
 *
 * \code
 * static
 * SCIP_RETCODE runShell(
 *   int                        argc,  
 *   char**                     argv,   
 *   const char*                defaultsetname 
 *   )
 * {
 *   SCIP* scip = NULL;
 *
 *   SCIP_CALL( SCIPcreate(&scip) );
 *
 *   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
 *
 *   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );
 *   SCIP_CALL( SCIPincludeEventHdlrBoundwriting(scip) );
 *
 *   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );
 *
 *   SCIP_CALL( SCIPfree(&scip) );
 *
 *   BMScheckEmptyMemory();
 *   
 *   return SCIP_OKAY;
 * }
 * \endcode
 *
 */
